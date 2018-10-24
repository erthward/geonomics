#!/usr/bin/python
#genome.py


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

#geonomics imports
from ops import mutation

#other imports
import numpy as np
import pandas as pd
from numpy import random as r
import matplotlib.pyplot as plt
from collections import OrderedDict as OD
import random
import bitarray


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

class _RecombinationPaths:
    def __init__(self, recomb_paths):
        self.recomb_paths = recomb_paths

    def _get_paths(self, n):
        return(random.sample(self.recomb_paths, n))

class Trait:
    def __init__(self, idx, name, phi, n_loci, mu, layer, alpha_distr_mu,
                 alpha_distr_sigma, gamma, univ_adv):
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
        self.gamma = gamma
        self.univ_adv = univ_adv

        self.loci = np.int64([])
        self.alpha = np.array([])


    def _get_phi(self, pop):
        if type(self.phi) in (float, int):
            phi = np.array([self.phi]*len(pop))
        else:
            phi = self.phi[pop.cells[:,1], pop.cells[:,0]]
        return(phi)

    def _set_loci(self, loci):
        self.loci = np.hstack((self.loci, np.array([*loci])))
        self.n_loci = self.loci.size

class GenomicArchitecture:
    def __init__(self, p, dom, r, g_params, land):
        #ploidy (NOTE: for now will be 2 by default; later could consider
        #enabling polyploidy)
        self.x = 2
        #total length (i.e. number of markers)
        self.L = g_params.L
        #length of each chromosome
        self.l_c = g_params.l_c
        #Dict of allele frequencies for the 1-alleles for all (numbered by
        #keys) chromosomes in haploid genome
        self.p = p
        #True/False regarding whether to allow a locus to affect the
        #phenotype of more than one trait; defaults to False
        self.pleiotropy  = g_params.pleiotropy
        #array of dominance values for all loci
        self.dom = dom
        #set the _use_dom attribute based on whether any loci have a 1 for
        #their dominance value
        self._use_dom = np.any(self.dom)
        self.sex = g_params.sex
        #Dict of recombination rates between each locus and the next, 
        #for all chroms 
        #(NOTE: first will be forced to 1/float(x), to effect independent
        #segregation of chroms, after which recomb occurs as a
        #crossing-over path down the chroms
        self.r = r
        self._n_recomb_paths_mem = g_params.n_recomb_paths_mem
        self._n_recomb_paths_tot = g_params.n_recomb_paths_tot

        #The recombination-paths object will be assigned here; used to
        #speed up large quantities of binomial draws needed for recombination
        self._recomb_paths = None
        #genome-wide neutral mutation rate
        self.mu_neut = g_params.mu_neut
        #a set to keep track of all loci that don't influence the 
        #phenotypes of any trait; defaults to all loci, then will be updated
        self.neut_loci = set(range(self.L))
        #a set to keep track of all loci that influence the phenotype of at
        #least one trait; after burn-in, will be updated
        self.nonneut_loci = set()

        #genome-wide deleterious mutation rate
        self.mu_delet = g_params.mu_delet
        self.delet_loci = OD()
        self.delet_alpha_distr_shape = g_params.delet_alpha_distr_shape
        self.delet_alpha_distr_scale = g_params.delet_alpha_distr_scale

        #add a dict of Trait objects, if necessary
        self.traits = None
        if 'traits' in [*g_params]:
            self.traits = _make_traits(g_params.traits, land)

        #a set containing eligible loci for mutation; after burn-in,
        #will be updated
        self._mutable_loci = set()
        #set self._mu_tot, the total per-site, per-generation mutation rate
        mus = [mu for mu in (self.mu_neut, self.mu_delet) if mu is not None]
        if self.traits is not None:
            mus = mus + [trt.mu for trt in self.traits.values()]
        self._mu_tot = sum(mus)
        self._mut_fns = self._make_mut_fns_dict()
        #set ._planned_muts to None, for now (this is not yet implemented,
        #but thinking about it
        self._planned_muts = None

    #method to make a _mut_fns dict, containing a function 
    #for each type of mutation for this population
    def _make_mut_fns_dict(self):
        mut_fns = {}
        if self.mu_neut > 0:
            fn = lambda pop,offspring: mutation._do_neutral_mutation(pop,
                                                                    offspring)
            mut_fns.update({'neut': fn})
        if self.mu_delet > 0:
            fn = lambda pop,offspring: mutation._do_deleterious_mutation(pop,
                                                                    offspring)
            mut_fns.update({'delet': fn})
        if self.traits is not None:
            for trait_num in self.traits:
                if self.traits[trait_num].mu > 0:
                    fn = lambda pop,offspring: mutation._do_trait_mutation(pop,
                                                        offspring, trait_num)
                    mut_fns.update({'t%i' % trait_num: fn})
        return mut_fns

    #method to draw mutation types for any number of mutations chosen
    #to occur in a given timestep 
    def _draw_mut_types(self, num):
        type_dict = {'neut': self.mu_neut,
                 'delet': self.mu_delet,
                **{'t%i' % (k):v.mu for k,v in self.traits.items()} }
        types = []
        probs = []
        for k,v in type_dict.items():
            types.append(k)
            probs.append(v)
        probs = [p/sum(probs) for p in probs]
        choices = r.choice(types, p = probs, size = num, replace = True)
        return(choices)

    #method for drawing an effect size for one or many loci 
    def _draw_trait_alpha(self, trait_num, n=1):
        alpha = r.normal(self.traits[trait_num].alpha_distr_mu,
                            self.traits[trait_num].alpha_distr_sigma, n)
        #set all effects to positive if the trait is monogenic
        #(because effects will be added to 0)
        if self.traits[trait_num].n_loci == 1:
            alpha = np.abs(alpha)
        return(alpha)

    #method for drawing new trait or deleterious loci (for a mutation)
    #from the currently neutral loci
    def _draw_mut_loci(self):
        loci = r.choice([*self._mutable_loci])
        #TODO: SHOULD I INCLUDE HERE A STATEMENT THAT CHECKS
        #IF len(gen_arch._mutable_loci) IS <= SOME SMALL
        #VALUE, AND IF SO FINDS FIXED SITES IN THE POPULATION AND ADDS
        #THEM TO IT??? (so that I'm not running that time-intensive
        #procedure often, but to make sure I'm almost certain not
        #to run out of mutable sites (unless I by chance draw more 
        #mutations in the next turn than there are remaining mutables))
        #DOES THIS IMPLY ANY STRANGE, POTENTIALLY PROBLEMATIC POPGEN ARTEFACTS?
        return(loci)

    #method for drawing new deleterious mutational fitness effects
    def _draw_delet_s(self):
        s = r.gamma(self.delet_alpha_distr_shape, self.delet_alpha_distr_scale)
        s = min(s, 1)
        return(s)

    #method for assigning loci to traits 
    def _set_trait_loci(self, trait_num, mutational=False,
                                            loci=None, alpha=None):
        #if this is not the result of a point mutation, but instead
        #either an initial setup or manually introduced, then grab the
        #number of loci to be assigned
        if mutational == False:
            n = self.traits[trait_num].n_loci
        #otherwise, assign a single locus
        else:
            n = 1
        if loci is not None:
            if not np.iterable(loci):
                loci = [loci]
            loci = set([*loci])
        #else, draw loci randomly, either allowing pleiotropy or not
        elif not self.pleiotropy:
            loci = set(r.choice([*self.neut_loci], size = n, replace = False))
        elif self.pleiotropy:
            loci = set(r.choice(range(self.L), size = n, replace = False))
        #update the trait's loci
        self.traits[trait_num]._set_loci(loci)
        #add these loci to self.non-neutral and remove from
        #self.neut_loci, to keep track of all loci underlying any
        #traits (for purpose of avoiding pleiotropy)
        self.nonneut_loci.update(loci)
        self.neut_loci.difference_update(loci)
        #if the effect size(s) is/are provided, use those
        if alpha is not None:
            if not np.iterable(alpha):
                alpha = np.array([alpha])
            effects = np.array([*alpha])
        #else, draw effects from a Gaussian dist with mean 0 and sigma
        #provided by trait params (as per e.g. Yeaman and Whitlock 2011)
        else:
            effects = self._draw_trait_alpha(trait_num, n)
        #check that loci and effects are of equal length
        assert len(loci) == len(effects), ('Lengths of the two arrays '
            'containing the new trait loci and their effects are not equal.')
        #then add the effects to the trait's alpha array
        self.traits[trait_num].alpha = np.hstack((self.traits[trait_num].alpha,
                                                                    effects))


    #method for creating and assigning the r_lookup attribute
    def _make_recomb_paths(self):
        self._recomb_paths = _RecombinationPaths(_make_recomb_paths_bitarrays(
                                                                        self))


    #method for plotting all allele frequencies for the population
    def _plot_allele_frequencies(self, pop):
        populome = np.stack([ind.genome for ind in pop.values()])
        freqs = populome.sum(axis = 2).sum(axis = 0)/(2*populome.shape[0])
        print(self.p.shape)
        print(freqs.shape)
        print(freqs)
        plt.plot(range(self.L), self.p, ':r')
        plt.plot(range(self.L), freqs, '-b')
        plt.show()


    #method for pickling a genomic architecture
    def write_pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#generate allele_freqs
def _draw_allele_freqs(l):
    return(r.beta(1,1,l))


#simulate genotypes
def _draw_genotype(p):
    return(r.binomial(1, p))


def _make_traits(traits_params, land):
    params_copy = {**traits_params}
    #get the number of traits to create
    num_traits = len(params_copy)
    #and set each Layer number using Layer names
    for k, v in params_copy.items():
        lyr_num = [num for num, lyr in land.items(
                                                ) if lyr.name == v.layer]
        assert len(lyr_num) == 1, ("Expected to find a single Layer with "
            "the Layer name indicated for Trait %s, but instead found "
            "%i.") % (k, len(lyr_num))
        v['layer'] = lyr_num[0]
    #then for each of i traits, unpack the ith components of the remaining
    #params to create the trait dict
    traits = {n: Trait(n, k_v[0], **k_v[1]) for n, k_v in enumerate(
                                                        params_copy.items())}
    return(traits)


#simulate linkage values
def _draw_r(g_params, recomb_rate_fn = None):

    #use custom recomb_fn, if provided
    if recomb_rate_fn != None:
        recomb_array = np.array([max(0,min(0.5,
                                    recomb_rate_fn())) for locus in range(L)])
        return(recomb_array)

    #otherwise, use default function with either default or custom param vals
    else:
        L = g_params.L

        param_vals = {'r_distr_alpha': 7e2, 'r_distr_beta': 7e3}

        for param in ['r_distr_alpha', 'r_distr_beta_']:
            if (param in g_params.keys() and g_params[param] is not None):
                param_vals[param] = g_params[param]

        recomb_array = np.array([max(0,min(0.5,
            recomb_rate)) for recomb_rate in r.beta(
            a = param_vals['r_distr_alpha'], b = param_vals['r_distr_beta'],
                                                            size = L)])
    #NOTE: for now, using the Beta, which can be very flexibly parameterized
    #NOTE: current default alpha/beta vals, after being subtracted from 0.5 
    #in sim_r function, will result in a tight distribution of r vals 
    #around 0.21 (~ 95% between 0.19 and 0.22)
    #NOTE: to essentially fix r at 0.5, Beta(1,1e7) will do...
    #NOTE: Ideally, would be good to look into to developing some sort
    #of mixture distribution to reasonably and flexibly model
    #map-distance values...

        return(recomb_array)


def _get_chrom_breakpoints(l_c, L):
    breakpoints = np.array([0]+list(np.cumsum(sorted(l_c))[:-1]))
    assert np.alltrue(np.diff(np.array(
        list(breakpoints) + [L])) == np.array(sorted(l_c))), ("The "
        "breakpoints assigned will not produce chromosomes of the "
        "correct length.")
    return(breakpoints)


#carry out recombination using the lookup array in a GenomicArchitecture object
def _make_recombinants(r_lookup, n_recombinants):
    recombinants = np.array([r.choice(r_lookup[i,],
        size = n_recombinants, replace = True) for i in range(len(r_lookup))])
    recombinants = np.cumsum(recombinants, axis = 0)%2
    return(recombinants)


def _make_recomb_array(g_params, recomb_values):
    #get L (num of loci) and l_c (if provided; num of loci per
    #chromsome) from the genome params dict
    L = g_params.L
    if ('l_c' in g_params.keys() and g_params['l_c'] != None
        and len(g_params['l_c']) > 1):
        l_c = g_params.l_c
        #and if l_c provided, check chrom lenghts sum to total number of loci
        assert sum(l_c) == L, ("The chromosome lengths provided do not sum to "
            "the number of loci provided.")
    else:
        l_c = [L]

    #if g_params.recomb_array (i.e a linkage map) manually provided (will
    #break if not a list, tuple, or np.array), then set that as the 
    #recomb_array, and check that len(recomb_array) == L
    if recomb_values is not None:
        assert len(recomb_values) == L, ("The length of the the table is the "
        "custom genomic-architecture file must be equal to the stipulated "
        "genome length ('L').")
        recomb_array = recomb_values[:]

    #otherwise, create recomb array
    else:
        #if a custom recomb_fn is provided, grab it
        if 'recomb_rate_custom_fn' in g_params.values():
            if g_params['recomb_rate_custom_fn'] is not None:
                recomb_rate_fn = g_params.recomb_rate_custom_fn
                assert callable(recomb_rate_fn), ("The 'recomb_rate_custom_fn'"
                    " provided in the parameters appears not to be defined "
                    "properly as a callable function.")
                #then call the _draw_r() function for each locus,
                #using custom recomb_fn
                recomb_array = _draw_r(g_params, recomb_fn = recomb_rate_fn)

        #otherwise, use the default _draw_r function to draw recomb rates
        else:
            recomb_array = _draw_r(g_params)

    #if more than one chromosome (i.e. if l_c provided in g_params dict and of
    #length >1), set recomb rate at the appropriate chrom breakpoints to 0.5
    if len(l_c) >1:
        bps = _get_chrom_breakpoints(l_c, L)
        recomb_array[bps] = 0.5
    #NOTE: Always set the first locus r = 0.5, to ensure independent 
    #assortment of homologous chromosomes
    recomb_array[0] = 0.5

    return(recomb_array, sorted(l_c))


#function to create a lookup array, for raster recombination of larger
#numbers of loci on the fly
#NOTE: size argument ultimately determines the minimum distance between
#probabilities (i.e. recombination rates) that can be modeled this way
def _make_recomb_paths_bitarrays(genomic_architecture,
                    lookup_array_size = 10000, n_recomb_paths_tot = 100000):

    if genomic_architecture._n_recomb_paths_mem is not None:
        lookup_array_size = genomic_architecture._n_recomb_paths_mem

    if genomic_architecture._n_recomb_paths_tot is not None:
        n_recomb_paths_tot = genomic_architecture._n_recomb_paths_tot

    lookup_array = np.zeros((len(genomic_architecture.r),lookup_array_size),
                                                            dtype = np.int8)

    for i, rate in enumerate(genomic_architecture.r):
        lookup_array[i,0:int(round(lookup_array_size*rate))] = 1

    recomb_paths = _make_recombinants(lookup_array, n_recomb_paths_tot).T
    bitarrays= tuple([_make_bitarray_recomb_subsetter(
                                                    p) for p in recomb_paths])

    return(bitarrays)


def _make_bitarray_recomb_subsetter(recomb_path):
    ba = bitarray.bitarray(list(recomb_path.reshape((recomb_path.size,))))
    ba_inv = bitarray.bitarray(list(np.array(ba) == False))
    tot = []
    for i in range(len(ba)):
        tot.append(ba[i])
        tot.append(ba_inv[i])
    return(bitarray.bitarray(tot))


#build the genomic architecture
def _make_genomic_architecture(pop_params, land):
    #get the genome parameters
    g_params = pop_params.gen_arch
    #get the custom genomic-architecture file, if provided
    gen_arch_file = None
    if 'gen_arch_file' in g_params.keys():
        if g_params.gen_arch_file is not None:
            gen_arch_file = pd.read_csv(g_params.gen_arch_file)
            assert len(gen_arch_file) == g_params.L, ('The custom '
            'genomic architecture file must contain a table with number of '
            'rows equal to the genome length stipulated in the genome '
            'parameters (params.comm[<pop_num>].gen_arch.L).')
            assert np.all(gen_arch_file['locus'].values ==
                          np.array([*range(len(gen_arch_file))])), ('The '
            "'locus' column of the custom genomic architecture file must "
            "contain serial integers from 0 to 1 - the length of the table.")
            assert (np.all(
              (gen_arch_file['dom'] == 0) + gen_arch_file['dom'] == 1)), ("The"
              " 'dom' column of the custom genomic architecture file must "
              "contain only 0s and 1s.")

    #also get the sex parameter and add it as an item in g_params
    g_params['sex'] = pop_params.mating.sex

    #draw locus-wise 1-allele frequencies, unless provided in 
    #custom gen-arch file
    if gen_arch_file is None:
        p = _draw_allele_freqs(g_params.L)
    else:
        p = gen_arch_file['p'].values

    #get the custom recomb_values, if provided in a custom gen-arch file
    recomb_values = None
    if gen_arch_file is not None:
        recomb_values = gen_arch_file['r'].values
        assert recomb_values[0] == 0.5, ("The top recombination rate value "
        "in the custom genomic-architecture file's 'r' column must be set "
        "to 0.5, to effect independent assortment of sister chromatids.")

    #set locus-wise dominance values for the 1-alleles, using the 'dom' value
    #in the gen_arch params, unless a gen_arch_file was provided
    if gen_arch_file is None:
        #create an L-length array of boolean integers (where 0 = codominance,
        #1 = dominance)
        dom = np.array([int(g_params.dom)] * g_params.L)
    else:
        #get the 'dom' column from the gen_arch_file
        dom = gen_arch_file['dom'].values

    r, l_c = _make_recomb_array(g_params, recomb_values = recomb_values)
    #in case g_params.l_c was missing or None, because only a single chromosome
    #is being simulated, then replace g_params.l_c with the returned l_c value
    g_params.l_c = l_c

    #now make the gen_arch object
    gen_arch = GenomicArchitecture(p, dom, r, g_params, land)

    #set the loci and effect sizes for each trait, using the custom gen-arch
    #file, if provided
    if gen_arch_file is not None:
        #convert the trait names in the 'trait' column of the file into 
        #their trait numbers (i.e. their keys in the gen_arch traits dict
        trait_names_nums = {
            trt.name: num for num, trt in gen_arch.traits.items()}
        gen_arch_file['trait'] = [trait_names_nums[
                val] for val in gen_arch_file['trait']]
        #turn the values in the 'trait' and 'alpha' columns into lists of
        #values, by splitting on commas
        #(this will allow people to assign a single locus
        #to more than one trait, i.e. to model pleiotropy
        gen_arch_file['trait'] = [
            [*map(int, row.split(','))] for row in gen_arch_file['trait']]
        gen_arch_file['alpha'] = [
            [*map(float, row.split(','))] for row in gen_arch_file['alpha']]
        #get the loci for each trait
        loci = {trt_num: np.array(range(len(gen_arch_file)))[
            [n in row for row in gen_arch_file[
            'trait']]] for n in gen_arch.traits.keys()} 
        #get the effect size for each locus
        alphas = {trt_num: np.concatenate([np.array(row[1]['alpha'])[
            [n == trt_num for n in row[1][
            'trait']]] for row in gen_arch_file.iterrows(
            )]) for trt_num in gen_arch.traits.keys()}
        #check that we got the same length of loci and effects for each trait
        for trt_num in loci.keys():
            assert len(loci[trt_num]) == len(alphas[trt_num]), ("Expected to "
                "receive the same number of loci and alphas (i.e. effect "
                "sizes) for trait number %i, but instead got %i loci and %i "
                "alphas.") % (trt_num, len(loci[trt_num]),
                len(alphas[trt_num]))
        #get the loci and effect sizes for each trait
        #trait_effects = {trait_num:
        #  {int(k): v for k,v in gen_arch_file[gen_arch_file['trait'] ==
        #                                trait_num][['locus','alpha']].values()}
        #    for trait_num in gen_arch.traits.keys()}
        #add the loci and effect sizes for each of the traits to the
        #Trait object in the GenomicArchitecture
        for trait_num in gen_arch.traits.keys():
            gen_arch._set_trait_loci(trait_num, mutational = False,
                                loci = loci[trait_num],
                                alpha = alphas[trait_num])
    #or else randomly set the loci and effect sizes for each trait
    else:
        if gen_arch.traits is not None:
            for trait_num in gen_arch.traits.keys():
                gen_arch._set_trait_loci(trait_num, mutational = False)

    assert len(set(range(gen_arch.L)).difference(
        gen_arch.neut_loci.union(gen_arch.nonneut_loci))) == 0, ("The union "
            "of the gen_arch.neut_loci and gen_arch.nonneut_loci sets does "
            "not contain all loci indicated by gen_arch.L")

    #create the r_lookup attribute
    gen_arch._make_recomb_paths()

    return gen_arch


#make genome
def _draw_genome(genomic_architecture):
    #NOTE: if for some reason any loci are not properly set to either 0 or 1,
    #they will stick out as 9s
    new_genome = np.ones([genomic_architecture.L, genomic_architecture.x],
                                                        dtype = np.int8) * 9
    for homologue in range(genomic_architecture.x):
        new_genome[:,homologue] = _draw_genotype(genomic_architecture.p)

    assert type(new_genome) == np.ndarray, ("A new genome must be an instance "
                                                            "of numpy.ndarray")
    assert np.shape(new_genome) == (genomic_architecture.L,
                            genomic_architecture.x), ("A new genome must "
                            "wind up with shape = (L, ploidy).")

    return(new_genome)


#function to reset genomes after burn-in
def _set_genomes(pop, burn_T, T):
    #use mean n_births at tail end of burn-in to estimate number of mutations,
    #and randomly choose set of neutral loci 
    #of that length to go into the pop.gen_arch._mutable_loci attribute
    n_muts = mutation._calc_estimated_total_mutations(pop, burn_T, T)
    #add a small number if n_muts evaluates to 0
    if n_muts == 0:
        neut_loci_remaining = pop.gen_arch.L - len(pop.gen_arch.nonneut_loci)
        n_muts = min(int(np.ceil(0.1 * neut_loci_remaining)),
                                                    neut_loci_remaining)
    muts = set(r.choice([*pop.gen_arch.neut_loci], n_muts, replace = False))
    pop.gen_arch._mutable_loci.update(muts)
    #set those loci's p values to 0 (i.e. non-segregating)
    pop.gen_arch.p[np.array([*muts])] = 0
    #now reassign genotypes to all individuals, using gen_arch.p
    [ind._set_genome(_draw_genome(pop.gen_arch)) for ind in pop.values()]
    #and then reset the individuals' phenotypes
    if pop.gen_arch.traits is not None:
        [ind._set_z(pop.gen_arch) for ind in pop.values()];


#method for loading a pickled genomic architecture
def read_pickled_genomic_architecture(filename):
    import cPickle
    with open(filename, 'rb') as f:
        gen_arch = cPickle.load(f)
    return gen_arch


