#!/usr/bin/python
# genome.py


'''
##########################################

Module name:              structs.genome

Module contents:          - definition of the Trait, _RecombinationPaths, and
                            GenomicArchitecture classes
                          - functions for simulation of genomic architecture,
                            simulation of a new genome, and other associated
                            functions


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               12-28-15
Documentation:            URL


##########################################
'''

# geonomics imports
from geonomics.ops import mutation

# other imports
import numpy as np
import pandas as pd
from numpy import random as r
import matplotlib.pyplot as plt
from collections import OrderedDict as OD
import warnings
import random
import bitarray

######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################


class _RecombinationPaths:
    def __init__(self, L, recomb_paths=None, fixed_r=None):
        self._L = L
        self.recomb_paths = recomb_paths
        self._fixed_r = fixed_r
        self._gen_paths_ad_hoc = self.recomb_paths is None

    def _get_paths(self, n):
        if self._gen_paths_ad_hoc:
            paths = _get_bitarray_subsetters_fixed_r(n, self._L, self._fixed_r)
        else:
            paths = random.sample(self.recomb_paths, n)
        return paths


class Trait:
    def __init__(self, idx, name, phi, n_loci, mu, layer, alpha_distr_mu,
                 alpha_distr_sigma, max_alpha_mag, gamma, univ_adv):
        self.idx = idx
        self.name = name
        self.phi = phi
        self.n_loci = n_loci
        if mu is None:
            mu = 0
        self.mu = mu
        self.lyr_num = layer
        self.alpha_distr_mu = alpha_distr_mu
        self.alpha_distr_sigma = alpha_distr_sigma
        self.max_alpha_mag = max_alpha_mag
        self.gamma = gamma
        self.univ_adv = univ_adv

        self.loci = np.int64([])
        self.alpha = np.array([])

    def _get_phi(self, spp):
        if type(self.phi) in (float, int):
            phi = np.array([self.phi]*len(spp))
        else:
            phi = self.phi[spp.cells[:, 1], spp.cells[:, 0]]
        return(phi)

    def _set_loci(self, loci):
        self.loci = np.hstack((self.loci, np.array([*loci])))
        self.n_loci = self.loci.size


class GenomicArchitecture:
    def __init__(self, p, dom, r, g_params, land):
        # ploidy (NOTE: for now will be 2 by default; later could consider
        # enabling polyploidy)
        self.x = 2
        # total length (i.e. number of markers)
        self.L = g_params.L
        # length of each chromosome
        self.l_c = g_params.l_c
        # Dict of allele frequencies for the 1-alleles for all (numbered by
        # keys) chromosomes in haploid genome
        self.p = p
        # True/False regarding whether to allow a locus to affect the
        # phenotype of more than one trait; defaults to False
        self.pleiotropy = g_params.pleiotropy
        # array of dominance values for all loci
        self.dom = dom
        # set the _use_dom attribute based on whether any loci have a 1 for
        # their dominance value
        self._use_dom = np.any(self.dom)
        self.sex = g_params.sex
        # Dict of recombination rates between each locus and the next,
        # for all chroms
        # (NOTE: first will be forced to 1/float(x), to effect independent
        # segregation of chroms, after which recomb occurs as a
        # crossing-over path down the chroms
        self.r = r
        self._n_recomb_paths_mem = g_params.n_recomb_paths_mem
        self._n_recomb_paths_tot = g_params.n_recomb_paths_tot

        # The recombination-paths object will be assigned here; used to
        # speed up large quantities of binomial draws needed for recombination
        self._recomb_paths = None
        # Get the allow_ad_hoc_recomb param, to determine whether or not to
        # allow the model to simulate recombination paths ad hoc (rather than
        # generate them at the beginning and then shuffle and draw from them
        # during the model run)
        # NOTE: this only works for homogeneous recombination across a genome
        # NOTE: this is significantly slower than pre-
        # generated recomb paths for combinations of large values of L (genome
        # size) and large values of N (mean pop size), but it allows combos of
        # very large values of those params to be run with less memory expense,
        # and it models exact recombination, rather than approximating it
        self._allow_ad_hoc_recomb = g_params.allow_ad_hoc_recomb
        # genome-wide neutral mutation rate
        self.mu_neut = g_params.mu_neut
        # a set to keep track of all loci that don't influence the
        # phenotypes of any trait; defaults to all loci, then will be updated
        self.neut_loci = set(range(self.L))
        # a set to keep track of all loci that influence the phenotype of at
        # least one trait; after burn-in, will be updated
        self.nonneut_loci = set()

        # genome-wide deleterious mutation rate
        self.mu_delet = g_params.mu_delet
        self.delet_loci = OD()
        self.delet_alpha_distr_shape = g_params.delet_alpha_distr_shape
        self.delet_alpha_distr_scale = g_params.delet_alpha_distr_scale

        # add a dict of Trait objects, if necessary
        self.traits = None
        if 'traits' in [*g_params]:
            self.traits = _make_traits(g_params.traits, land)

        # a set containing eligible loci for mutation; after burn-in,
        # will be updated
        self._mutable_loci = set()
        # set self._mu_tot, the total per-site, per-generation mutation rate
        mus = [mu for mu in (self.mu_neut, self.mu_delet) if mu is not None]
        if self.traits is not None:
            mus = mus + [trt.mu for trt in self.traits.values()]
        self._mu_tot = sum(mus)
        # attribute that will be replaced with the estimated total number of
        # mutations per iteration, once it's estimated at the end of the first
        # burn-in
        self._est_tot_muts_per_it = None
        self._mut_fns = self._make_mut_fns_dict()
        # set ._planned_muts to None, for now (this is not yet implemented,
        # but thinking about it
        self._planned_muts = None

    # method to make a _mut_fns dict, containing a function
    # for each type of mutation for this species
    def _make_mut_fns_dict(self):
        mut_fns = {}
        if self.mu_neut > 0:
            def fn(spp, offspring):
                mutation._do_neutral_mutation(spp, offspring)
            mut_fns.update({'neut': fn})
        if self.mu_delet > 0:
            def fn(spp, offspring):
                mutation._do_deleterious_mutation(spp, offspring)
            mut_fns.update({'delet': fn})
        if self.traits is not None:
            for trait_num in self.traits:
                if self.traits[trait_num].mu > 0:
                    def fn(spp, offspring, trait_num=trait_num):
                        mutation._do_trait_mutation(spp, offspring, trait_num)
                    mut_fns.update({'t%i' % trait_num: fn})
        return mut_fns

    # method to draw mutation types for any number of mutations chosen
    # to occur in a given timestep
    def _draw_mut_types(self, num):
        type_dict = {'neut': self.mu_neut,
                     'delet': self.mu_delet}
        if self.traits is not None:
            [type_dict.update(item) for item in {'t%i' % (
                            k): v.mu for k, v in self.traits.items()}.items()]
        types = []
        probs = []
        for k, v in type_dict.items():
            types.append(k)
            probs.append(v)
        probs = [p/sum(probs) for p in probs]
        choices = r.choice(types, p=probs, size=num, replace=True)
        return(choices)

    # method for drawing an effect size for one or many loci
    def _draw_trait_alpha(self, trait_num, n=1):
        mu = self.traits[trait_num].alpha_distr_mu
        sigma = self.traits[trait_num].alpha_distr_sigma
        max_alpha_mag = self.traits[trait_num].max_alpha_mag
        if max_alpha_mag is not None:
            min_alpha = -1 * max_alpha_mag
            max_alpha = max_alpha_mag
        else:
            min_alpha = max_alpha = max_alpha_mag
        # use mu value as the fixed effect size, if sigma is 0
        if sigma == 0:
            alpha = mu * np.array([1 - (i % 2)*2 for i in range(n)])
        else:
            alpha = r.normal(self.traits[trait_num].alpha_distr_mu,
                             self.traits[trait_num].alpha_distr_sigma, n)
            if max_alpha_mag is not None:
                alpha = np.clip(alpha, min_alpha, max_alpha)
        # otherwise use mu and sigma to draw effects from
        # (because effects will be added to 0)
        if self.traits[trait_num].n_loci == 1:
            alpha = np.abs(alpha)
        return(alpha)

    # method for drawing new deleterious mutational fitness effects
    def _draw_delet_s(self):
        s = r.gamma(self.delet_alpha_distr_shape, self.delet_alpha_distr_scale)
        s = min(s, 1)
        return(s)

    # method for assigning loci to traits
    def _set_trait_loci(self, trait_num, mutational=False,
                        loci=None, alpha=None):
        # if this is not the result of a point mutation, but instead
        # either an initial setup or manually introduced, then grab the
        # number of loci to be assigned
        if mutational is False:
            n = self.traits[trait_num].n_loci
            assert n <= self.L, ("The number of loci parameterized for "
                                 "trait number %i ('n_loci') is greater "
                                 "than the length of the genome!")
        # otherwise, assign a single locus
        else:
            n = 1
        if loci is not None:
            if not np.iterable(loci):
                loci = [loci]
            assert len(set([*loci])) == len(loci), ("Some of the trait "
                                                    "loci provided appear "
                                                    "to be repeated.")
        # else, draw loci randomly, either allowing pleiotropy or not
        elif not self.pleiotropy:
            loci = set(r.choice([*self.neut_loci], size=n, replace=False))
        elif self.pleiotropy:
            loci = set(r.choice(range(self.L), size=n, replace=False))
        # update the trait's loci
        self.traits[trait_num]._set_loci(loci)
        # add these loci to self.non-neutral and remove from
        # self.neut_loci, to keep track of all loci underlying any
        # traits (for purpose of avoiding pleiotropy)
        self.nonneut_loci.update(loci)
        self.neut_loci.difference_update(loci)
        # if the effect size(s) is/are provided, use those
        if alpha is not None:
            if not np.iterable(alpha):
                alpha = np.array([alpha])
            effects = np.array([*alpha])
        # else, draw effects from a Gaussian dist with mean 0 and sigma
        # provided by trait params (as per e.g. Yeaman and Whitlock 2011)
        else:
            effects = self._draw_trait_alpha(trait_num, n)
        # if this is drawing the effect of a monogenic trait's locus (i.e. if
        # it's drawing only one value, but it's not mutational), then coerce
        # the effects to [0.5]
        if not mutational and n == 1:
            effects = np.array([0.5])
        # check that loci and effects are of equal length
        assert len(loci) == len(effects), ('Lengths of the two arrays '
                                           'containing the new trait loci and '
                                           'their effects are not equal.')
        # then add the effects to the trait's alpha array
        self.traits[trait_num].alpha = np.hstack((self.traits[trait_num].alpha,
                                                  effects))

    # method for creating and assigning the r_lookup attribute
    def _make_recomb_paths(self):
        self._recomb_paths = _RecombinationPaths(self.L,
                                                 *_make_recomb_paths_bitarrays(
                                                                        self))

    # method for plotting all allele frequencies for the species
    def _plot_allele_frequencies(self, spp):
        speciome = np.stack([ind.g for ind in spp.values()])
        freqs = speciome.sum(axis=2).sum(axis=0) / (2*speciome.shape[0])
        plt.plot(range(self.L), self.p, ':r')
        plt.plot(range(self.L), freqs, '-b')
        plt.show()

    # method for pickling a genomic architecture
    def write_pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

# generate allele_freqs
def _draw_allele_freqs(l):
    return(r.beta(1, 1, l))


# simulate genotypes
def _draw_genotype(p):
    return(r.binomial(1, p))


def _make_traits(traits_params, land):
    params_copy = {**traits_params}
    # and set each Layer number using Layer names
    for k, v in params_copy.items():
        # the first time this is run during a model with random
        # communities at each iteration, the layer identified in
        # the traits_params will be a string indicating the layer's name
        if isinstance(v.layer, str):
            lyr_num = [num for num, lyr in land.items(
                                                ) if lyr.name == v.layer]
        # the second and later times it's run during such a model (see
        # previous comment), it will already have been swapped out for an int
        # indicating the layer's index num
        elif isinstance(v.layer, int):
            lyr_num = [num for num, lyr in land.items(
                                                ) if lyr.idx == v.layer]
        assert len(lyr_num) == 1, ("Expected to find a single Layer with "
                                   "the Layer name indicated for Trait %s, "
                                   "but instead found "
                                   "%i.") % (k, len(lyr_num))
        v['layer'] = lyr_num[0]
    # then for each of i traits, unpack the ith components of the remaining
    # params to create the trait dict
    traits = {n: Trait(n, k_v[0], **k_v[1]) for n, k_v in enumerate(
                                                        params_copy.items())}
    # for all monogenic traits, if the trait doesn't already have a 0 mutation
    # rate then coerce it to 0
    # NOTE: this is the because there is a good reason to use 0 as the 
    # baseline phenotype for monogenic traits but 0.5 as the baseline
    # for multilocus traits (to preserve symmetry of ability for individuals
    # to have phenotypes beyond 0 and 1 for a multilocus trait and thus
    # experience stabilizing selection within their optimal habitat), but
    # this means that there would be a huge problem if a monogenic trait
    # underwent an adaptive mutation and became polyenic bceause all of the
    # individuals' phenotypes would suddenly have to be recalculated and would
    # suddenly completely change
    for n, trt in traits.items():
        if trt.n_loci == 1:
            if trt.mu != 0:
                warnings.warn(("Coercing Trait %i ('%s') to a "
                               "0 mutation rate because it is monogenic."))
                trt.mu = 0
    return(traits)


# simulate linkage values
def _draw_r(g_params, recomb_rate_fn=None):

    # use custom recomb_fn, if provided
    if recomb_rate_fn is not None:
        recomb_array = np.array([max(0, min(0.5,
                                recomb_rate_fn())) for locus in range(
                                g_params.L)])
        return(recomb_array)

    # otherwise, use default function with either default or custom param vals
    else:
        L = g_params.L

        # if either or both distribution parameters are None, then set all
        # recomb rates to 0.5
        if (g_params.r_distr_alpha is None and g_params.r_distr_beta is None):
            recomb_array = np.array([0.5]*L)
        # or fix the recomb rates at r_distr_alpha if it has a numeric value
        # but r_distr_beta is None
        elif ((g_params.r_distr_alpha is not None)
              and (g_params.r_distr_beta is None)):
            recomb_array = np.array([g_params.r_distr_alpha]*L)
        # or throw an error if r_distr_beta is non-None but alpha is None
        elif ((g_params.r_distr_alpha is None)
              and (g_params.r_distr_beta is not None)):
            raise ValueError(("The genomic architecture's 'r_distr_beta' "
                              "argument cannot have a numeric value if "
                              "'r_distr_alpha' is None."))
        # or if they both have numeric values, then draw the recombination
        # rates from a beta distribution
        else:
            recomb_array = np.clip(r.beta(a=g_params.r_distr_alpha,
                                   b=g_params.r_distr_beta, size=L),
                                   a_min=0, a_max=0.5)

        return(recomb_array)


def _get_chrom_breakpoints(l_c, L):
    breakpoints = np.array([0]+list(np.cumsum(sorted(l_c))[:-1]))
    assert_msg = ("The breakpoints assigned will not produce chromosomes "
                  "of the correct length.")
    assert np.alltrue(np.diff(np.array(
        list(breakpoints) + [L])) == np.array(sorted(l_c))), assert_msg
    return(breakpoints)


# carry out recombination using lookup array in a GenomicArchitecture object
def _make_recombinants(r_lookup, n_recombinants):
    recombinants = np.array([r.choice(r_lookup[i, ],
                            size=n_recombinants,
                            replace=True) for i in range(len(r_lookup))])
    recombinants = np.cumsum(recombinants, axis=0) % 2
    return(recombinants)


# method to simulate recombination for a genomic arch with a fixed recomb rate
def _get_bitarray_subsetters_fixed_r(n, L, fixed_r):
    # TODO: make this work for n individs at once!
    indep_assort = r.binomial(n=1, p=0.5, size=n).reshape((n, 1))
    n_recombs = r.binomial(n=n*(L-1), p=fixed_r)
    sites = r.choice(a=range(n*(L-1)), size=n_recombs, replace=False)
    recombs = np.zeros(n*(L-1))
    recombs[sites] = 1
    recombs = np.cumsum(recombs.reshape((n,L-1)), axis=1) % 2
    paths = np.hstack((indep_assort, (recombs + indep_assort) % 2))
    subsetters = [_make_bitarray_recomb_subsetter(path) for path in paths]
    #conglom_subsetter = _make_bitarray_recomb_subsetter(paths.flatten())
    #subsetters = [*np.array(conglom_subsetter).reshape((n, 2 * L))]
    return subsetters


def _make_recomb_array(g_params, recomb_values):
    # get L (num of loci) and l_c (if provided; num of loci per
    # chromsome) from the genome params dict
    L = g_params.L
    if ('l_c' in g_params.keys() and g_params['l_c'] is None
            and len(g_params['l_c']) > 1):
        l_c = g_params.l_c
        # and if l_c provided, check chrom lenghts sum to total number of loci
        assert sum(l_c) == L, ("The chromosome lengths provided do not sum to "
                               "the number of loci provided.")
    else:
        l_c = [L]

    # if g_params.recomb_array (i.e a linkage map) manually provided (will
    # break if not a list, tuple, or np.array), then set that as the
    # recomb_array, and check that len(recomb_array) == L
    if recomb_values is not None:
        assert len(recomb_values) == L, ("The length of the the table is the "
                                         "custom genomic-architecture file "
                                         "must be equal to the stipulated "
                                         "genome length ('L').")
        recomb_array = recomb_values[:]

    # otherwise, create recomb array
    else:
        # if a custom recomb_fn is provided, grab it
        if 'recomb_rate_custom_fn' in g_params.values():
            if g_params['recomb_rate_custom_fn'] is not None:
                recomb_rate_fn = g_params.recomb_rate_custom_fn
                assert callable(recomb_rate_fn), ("The 'recomb_rate_custom_fn'"
                                                  " provided in the "
                                                  "parameters appears not "
                                                  "to be defined properly as "
                                                  "a callable function.")
                # then call the _draw_r() function for each locus,
                # using custom recomb_fn
                recomb_array = _draw_r(g_params, recomb_fn=recomb_rate_fn)

        # otherwise, use the default _draw_r function to draw recomb rates
        else:
            recomb_array = _draw_r(g_params)

    # if more than one chromosome (i.e. if l_c provided in g_params dict and of
    # length >1), set recomb rate at the appropriate chrom breakpoints to 0.5
    if len(l_c) > 1:
        bps = _get_chrom_breakpoints(l_c, L)
        recomb_array[bps] = 0.5
    # NOTE: Always set the first locus r = 0.5, to ensure independent
    # assortment of homologous chromosomes
    recomb_array[0] = 0.5

    return(recomb_array, sorted(l_c))


# function to create a lookup array, for raster recombination of larger
# numbers of loci on the fly
# NOTE: size argument ultimately determines the minimum distance between
# probabilities (i.e. recombination rates) that can be modeled this way
def _make_recomb_paths_bitarrays(genomic_architecture,
                                 lookup_array_size=10000,
                                 n_recomb_paths_tot=100000):
    # only make bitarrays if recombination rates are heterogeneous
    if (len(np.unique(genomic_architecture.r[1:])) > 1
        or not genomic_architecture._allow_ad_hoc_recomb):

        if genomic_architecture._n_recomb_paths_mem is not None:
            lookup_array_size = genomic_architecture._n_recomb_paths_mem

        if genomic_architecture._n_recomb_paths_tot is not None:
            n_recomb_paths_tot = genomic_architecture._n_recomb_paths_tot

        # raise a warning if the smallest recombination rate is smaller
        # than the precision that can be modeled using the chosen number
        # of recombination paths to be held in memory (i.e. n_recomb_paths_mem,
        # in the params file; lookup_array_size in the arguments here)
        if 1/lookup_array_size > min(genomic_architecture.r):
            warnings.warn(("The number of recombination paths to be held in "
                           "memory (i.e. parameter 'n_recomb_paths_mem' in the"
                           " parameters file) provides for less precision in "
                           "Geonomics' approximation of recombination rates "
                           "than the minimum non-zero recombination rate "
                           "stipulated in your genomic architecture. (The "
                           "precision of this estimation is determined by "
                           "1/n_recomb_paths_mem.) Consider either increasing "
                           "n_recomb_paths_mem or increasing your minimum "
                           "non-zero recombination rate."))

        lookup_array = np.zeros((len(genomic_architecture.r),
                                 lookup_array_size), dtype=np.int8)

        for i, rate in enumerate(genomic_architecture.r):
            # NOTE: taking the max of the count within lookup_array_size that
            # represents the recombination rate and the integer of testing
            # rate != 0 ensures that for every nonzero recombination rate
            # we get at least a single recombination path that recombines
            # at that locus
            lookup_array[i, 0:max(int(round(lookup_array_size * rate)),
                                  int(rate != 0))] = 1

        recomb_paths = _make_recombinants(lookup_array, n_recomb_paths_tot).T
        bitarrays = tuple([_make_bitarray_recomb_subsetter(
                                                    p) for p in recomb_paths])

        # and set fixed_r to None
        fixed_r = None

    # if recombination rates are homoegenous, then just return None for the
    # bitarrays, and return the fixed recombination rate as fixed_r, because
    # recombinants will be quickly generated on the fly
    else:
        fixed_r = np.unique(genomic_architecture.r)[0]
        bitarrays = None

    return bitarrays, fixed_r


def _make_bitarray_recomb_subsetter(recomb_path):
    chrom1 = bitarray.bitarray([*recomb_path.reshape((recomb_path.size,))])
    chrom0 = chrom1[:]
    chrom0.invert()
    #chrom0 = bitarray.bitarray([*np.invert(chrom1)])
    # NOTE: This will create a binary subsetter than can be used to subset a
    # 2 x L genome that's been flatted into a 1 x 2L vector (hence the
    # intercalated use of ba and ba_inv)
    subsetter = bitarray.bitarray([val for item in [*zip(
                                            chrom0, chrom1)] for val in item])
    return(subsetter)


# build the genomic architecture
def _make_genomic_architecture(spp_params, land):
    # get the genome parameters
    g_params = spp_params.gen_arch
    # get the custom genomic-architecture file, if provided
    gen_arch_file = None
    if 'gen_arch_file' in g_params.keys():
        if g_params.gen_arch_file is not None:
            gen_arch_file = pd.read_csv(g_params.gen_arch_file)
            # ensure that trait and alpha columns are strings to start (because
            # in many cases alpha would likely read in as a float, if each row
            # has at most one alpha value because there's no pleiotropy, but in
            # case of pleiotropy I need to be able to use the str.split method
            # on each row's values further down)
            gen_arch_file['trait'] = [str(v) for v in gen_arch_file['trait']]
            gen_arch_file['alpha'] = [str(v) for v in gen_arch_file['alpha']]
            assert_msg = ('The custom genomic architecture file must contain '
                          'a table with number of rows equal to the genome '
                          'length stipulated in the genome parameters '
                          '(params.comm[<spp_name>].gen_arch.L).')

            assert len(gen_arch_file) == g_params.L, assert_msg
            assert_msg = ("The 'locus' column of the custom genomic "
                          "architecture file must contain serial integers "
                          "from 0 to 1 - the length of the table.")
            assert np.all(gen_arch_file['locus'].values ==
                          np.array([*range(len(gen_arch_file))])), assert_msg
            assert (np.all(
                (gen_arch_file['dom'] == 0) + gen_arch_file['dom'] == 1)), (
                "The 'dom' column of the custom genomic architecture file "
                "must contain only 0s and 1s (where 0 indicates codominance "
                "of the 0 and 1 alleles, 1 indicates that the 1 allele "
                "is dominant).")
            # check that each trait has in that file the number of loci
            # indicated by n_loci in that trait's params dict
            if 'traits' in [*g_params]:
                all_traits = [trt.strip() for row in [val.split(
                    ',') for val in gen_arch_file['trait']] for trt in row]
                for trt_name, trt in g_params.traits.items():
                    n_loci_in_file = sum(
                                    [trt == trt_name for trt in all_traits])
                    assert_msg = ("The number of times a Trait is appears "
                                  "in the custom genomic architecture file "
                                  "must be equivalent to the number of loci "
                                  "subtending that Trait as indicated by the "
                                  "'n_loci' key in its section of the "
                                  "parameters file.")

                    assert n_loci_in_file == trt.n_loci, assert_msg

    # also get the sex parameter and add it as an item in g_params
    g_params['sex'] = spp_params.mating.sex

    # draw locus-wise 1-allele frequencies, unless provided in
    # custom gen-arch file
    if gen_arch_file is None:
        if g_params.start_p_fixed:
            p = np.array([0.5]*g_params.L)
        else:
            p = _draw_allele_freqs(g_params.L)
    else:
        p = gen_arch_file['p'].values

    # get the custom recomb_values, if provided in a custom gen-arch file
    recomb_values = None
    if gen_arch_file is not None:
        recomb_values = gen_arch_file['r'].values
        assert recomb_values[0] == 0.5, ("The top recombination rate value "
                                         "in the custom genomic-architecture "
                                         "file's 'r' column must be set "
                                         "to 0.5, to effect independent "
                                         "assortment of sister chromatids.")

    # set locus-wise dominance values for the 1-alleles, using the 'dom' value
    # in the gen_arch params, unless a gen_arch_file was provided
    if gen_arch_file is None:
        # create an L-length array of boolean integers (where 0 = codominance,
        # 1 = dominance)
        dom = np.array([int(g_params.dom)] * g_params.L)
    else:
        # get the 'dom' column from the gen_arch_file
        dom = gen_arch_file['dom'].values

    r, l_c = _make_recomb_array(g_params, recomb_values=recomb_values)
    # in case g_params.l_c was missing or None, because only a single
    # chromosome is being simulated, then replace g_params.l_c with
    # the returned l_c value
    g_params.l_c = l_c

    # now make the gen_arch object
    gen_arch = GenomicArchitecture(p, dom, r, g_params, land)

    # set the loci and effect sizes for each trait, using the custom gen-arch
    # file, if provided
    if gen_arch_file is not None:
        # convert the trait names in the 'trait' column of the file into
        # lists of their trait numbers (i.e. their keys in the gen_arch
        # traits dict)
        trait_names_nums = {
            trt.name: num for num, trt in gen_arch.traits.items()}
        gen_arch_file['trait'] = [[trait_names_nums[
            val] for val in [x.strip() for x in row.split(
                                    ',')]] for row in gen_arch_file['trait']]
        # turn the values in the 'alpha' column into lists of
        # values, by splitting on commas
        # (this will allow people to assign a single locus
        # to more than one trait, i.e. to model pleiotropy)
        gen_arch_file['alpha'] = [
            [*map(float, row.split(','))] for row in gen_arch_file['alpha']]
        # get the loci and effect sizes for each trait
        loci = {}
        alphas = {}
        for trt_num in gen_arch.traits.keys():
            loci[trt_num] = np.array([*gen_arch_file['locus']])[
                [trt_num in row for row in gen_arch_file['trait']]]
            alphas[trt_num] = np.concatenate([np.array(row['alpha'])[
                [n == trt_num for n in row[
                    'trait']]] for i, row in gen_arch_file.iterrows()])
        # check that we got the same length of loci and effects for each trait
        for trt_num in loci.keys():
            assert_msg = ("Expected to receive the same number of loci and "
                          "alphas (i.e. effect sizes) for trait number %i, "
                          "but instead got %i loci and %i "
                          "alphas.") % (trt_num, len(loci[trt_num]),
                                        len(alphas[trt_num]))

            assert len(loci[trt_num]) == len(alphas[trt_num]), assert_msg

        # add the loci and effect sizes for each of the traits to the
        # Trait object in the GenomicArchitecture
        for trait_num in gen_arch.traits.keys():
            gen_arch._set_trait_loci(trait_num, mutational=False,
                                     loci=loci[trait_num],
                                     alpha=alphas[trait_num])
    # or else randomly set the loci and effect sizes for each trait
    else:
        if gen_arch.traits is not None:
            for trait_num in gen_arch.traits.keys():
                gen_arch._set_trait_loci(trait_num, mutational=False)

    assert_msg = ("The union of the gen_arch.neut_loci and "
                  "gen_arch.nonneut_loci sets does not contain all loci "
                  "indicated by gen_arch.L")
    assert len(set(range(gen_arch.L)).difference(
        gen_arch.neut_loci.union(gen_arch.nonneut_loci))) == 0, assert_msg

    # create the r_lookup attribute
    gen_arch._make_recomb_paths()

    return gen_arch


# make genome
def _draw_genome(genomic_architecture):
    # NOTE: if for some reason any loci are not properly set to either 0 or 1,
    # they will stick out as 9s
    new_genome = np.ones([genomic_architecture.L, genomic_architecture.x],
                         dtype=np.int8) * 9
    for homologue in range(genomic_architecture.x):
        new_genome[:, homologue] = _draw_genotype(genomic_architecture.p)

    assert type(new_genome) == np.ndarray, ("A new genome must be an instance "
                                            "of numpy.ndarray")
    assert np.shape(new_genome) == (genomic_architecture.L,
                                    genomic_architecture.x), ("A new genome "
                                                              "must wind up "
                                                              "with shape = "
                                                              "(L, ploidy).")

    return(new_genome)


# check whether a GenomicArchitecture has enough mutable loci to provide
# for all the mutations expected during a Model run (i.e. iteration)
def _check_enough_mutable_loci(spp, burn_T, T):
    # grab the estimated total number of muts per iteration, if it's already
    # been estimated
    if spp.gen_arch._est_tot_muts_per_it is not None:
        n_muts = spp.gen_arch._est_tot_muts_per_it
    # otherwise estimate it, then store it at spp.gen_arch._est_tot_muts_per_it
    else:
        n_muts = mutation._calc_estimated_total_mutations(spp, burn_T, T)
        spp.gen_arch._est_tot_muts_per_it = n_muts
    neut_loci_remaining = spp.gen_arch.L - len(spp.gen_arch.nonneut_loci)
    # throw error if n_muts is greater than the number of
    # available mutable loci
    if n_muts > neut_loci_remaining:
        raise ValueError(("Before a Model is run, Geonomics sets aside a set "
                          "of loci where mutations can occur (known as "
                          "'mutable loci'). The mutation rates in this "
                          "parameterization imply that at least %i mutable "
                          "loci will be needed to provide space for all "
                          "mutations expected over the course of the each "
                          "Model iteration, but only %i non-neutral loci are "
                          "available. Either the mutation rates must be "
                          "reduced, the runtime of each iteration must be "
                          "reduced, the genome length must be increased, or "
                          "some combination thereof.") % (n_muts,
                                                          neut_loci_remaining))


# function to reset genomes after burn-in
def _set_genomes(spp, burn_T, T):
    # use mean n_births at tail end of burn-in to estimate number of mutations,
    # and randomly choose set of neutral loci
    # of that length to go into the spp.gen_arch._mutable_loci attribute
    n_muts = mutation._calc_estimated_total_mutations(spp, burn_T, T)
    # add a small number if n_muts evaluates to 0
    neut_loci_remaining = spp.gen_arch.L - len(spp.gen_arch.nonneut_loci)
    if n_muts == 0:
        n_muts = min(int(np.ceil(0.1 * neut_loci_remaining)),
                     neut_loci_remaining)
    n_muts = min(n_muts, neut_loci_remaining)

    # skip this step and force the neutral mutation rate to 0, if there are no
    # neutral loci in the genome as it was configured
    if len(spp.gen_arch.neut_loci) == 0:
        warnings.warn(('This species has been parameterized without any'
               ' neutral loci, leaving no place for mutations (neutral or not)'
               ' to land. Thus all mutation rates will be forced to 0.'))
        spp.gen_arch.mu_neut = 0
        spp.gen_arch.mu_delet = 0
        for trt in spp.gen_arch.traits.values():
            trt.mu = 0
    # or just outright skip this step if the parameterization included only
    # mutation rates of 0
    elif spp.gen_arch._mu_tot == 0:
        pass
    # otherwise choose and update mutable loci
    else:
        muts = set(r.choice([*spp.gen_arch.neut_loci], n_muts, replace=False))
        spp.gen_arch._mutable_loci.update(muts)
        # set those loci's p values to 0 (i.e. non-segregating)
        spp.gen_arch.p[np.array([*muts])] = 0
    # now reassign genotypes to all individuals, using gen_arch.p
    [ind._set_g(_draw_genome(spp.gen_arch)) for ind in spp.values()]
    # and then reset the individuals' phenotypes
    if spp.gen_arch.traits is not None:
        [ind._set_z(spp.gen_arch) for ind in spp.values()]


# method for loading a pickled genomic architecture
def read_pickled_genomic_architecture(filename):
    import cPickle
    with open(filename, 'rb') as f:
        gen_arch = cPickle.load(f)
    return gen_arch
