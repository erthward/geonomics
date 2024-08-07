#!/usr/bin/python
# genome.py

# flake8: noqa


'''
Classes, associated methods, and supporting functions for all genomic
components
'''

# geonomics imports
from geonomics.ops import mutation
from geonomics.structs.individual import Individual
from geonomics.utils.viz import _check_display

# other imports
import numpy as np
import pandas as pd
from numpy import random as r
import matplotlib as mpl
_check_display()
import matplotlib.pyplot as plt
from collections import OrderedDict as OD
import warnings
import random
import bisect
import bitarray
import msprime
import tskit


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################


# an error to raise if a simulation is not given enough neutral loci to fit
# the total number of expected mutations over the total parameterized runtime
# under the infinite sites model
class MutationRateError(Exception):
    pass


class Recombinations:
    """
    Container for the parameters and data necessary to simulate recombination.

    Not intended for public use.
    """
    def __init__(self, L, positions, n, r_distr_alpha, r_distr_beta,
                 recomb_rates, jitter_breakpoints):
        # define the bitarray unit corresponding to each homologue
        self._bitarray_dict= {0: '10', 1: '01'}

        # genome length
        self._L = L
        # organize the potential recombination breakpoint positions
        if positions is None:
            positions = np.arange(self._L)
        else:
            positions.sort()
        positions = [*positions]
        # optimize the data-type, to save some memory
        if len(positions) <= 2**16:
            self._positions = np.int16(np.array(positions))
        else:
            self._positions = np.int32(np.array(positions))

        # number of recombination events to simulate and cache for the model
        # NOTE: the higher this number, the more accurately geonomics
        #       will simulate the stipulated recombination rates
        self._n = n
        # alpha and beta parameters for the beta distribution from which
        # recombination rates can be drawn
        self._r_distr_alpha = r_distr_alpha
        self._r_distr_beta = r_distr_beta
        # flag indicating whether or not breakpoint locations should be
        # jittered slightly off of their x.5 defaults, so that
        # tskit.TreeSequence will correctly track and report the number of trees
        self._jitter_breakpoints = jitter_breakpoints
        # get the recombination rates
        if recomb_rates is not None:
            assert len(recomb_rates) == len(positions), ("Lengths of provided "
                                                         "recombination "
                                                         "rates and "
                                                         "recombination "
                                                         "positions don't "
                                                         "match!")
            assert recomb_rates[0] == 0, ("The first recombination rate (i.e. "
                                          "index 0) must be 0 (because no "
                                          "recombination can be modeled as "
                                          "taking place prior to the first "
                                          "locus being simulated).")
            self._rates = recomb_rates
        else:
            self._rates = self._draw_recombination_rates()

    def _set_events(self, use_subsetters, nonneutral_loci, use_tskit):
        # get the potential recombination breakpoints, their recombination
        # rates, and the cache of simulated recombination events to be used
        # by the model
        # NOTE: each item will only be generated and stored if it is needed,
        #       depending on whether or not tskit is being used and, if it is
        #       being used, whether or not there are non-neutral loci whose
        #       genotypes need to be tracked
        breakpoints, subsetters = self._draw_recombination_events(
                                            use_subsetters=use_subsetters,
                                            nonneutral_loci=nonneutral_loci,
                                            use_tskit=use_tskit)
        # if using tskit, set seg_info and use the breakpoints dict
        if use_tskit:
            self._breakpoints = breakpoints
            self._set_seg_info()
        # otherwise don't set seg_info and make the breakpoints dict None
        else:
            self._breakpoints = None
        # either way, set the subsetters
        self._subsetters = subsetters

    def _get_events(self, size):
        events = random.sample(self._events, size)
        return events

    def _get_subsetter(self, event_key):
        return self._subsetters[event_key]


    # update all the recombination subsetters in accord with a new non-neutral
    # mutation
    def _update_subsetters(self, mutation_loc, genotype_arr_mutation_idx):
        for k in self._subsetters.keys():
            # determine whether the number of crossovers that has happened up
            # to the mutation's locus has been even (in which case the
            # recomb path is 'back' on homologue 0) or odd (in which
            # case it is on homologue 1)
            # NOTE: this is because all recomb paths start on homologue 0,
            #       then are used either for subsetting starting from either
            #       homol 0 or homol 1 on the fly in mating.py using the 
            #       variable start_homologues
            #       in function _do_mating_singl_offspring
            subsetter_homologue = bisect.bisect_left(self._breakpoints[k],
                                                     mutation_loc) % 2
            # use the bitarray_dict to get the 2-digit binary string
            # corresponding to the homologue that the recomb path is on
            # when it reaches mutation_loc, to be inserted into the subsetter
            # in the appropriate location
            subsetter_insert = self._bitarray_dict[subsetter_homologue]
            # find the insertion-point index of the current subsetter 
            # (just twice the non-neutral genotype array's index)
            subsetter_idx = genotype_arr_mutation_idx * 2
            # create the new subsetter by concatenating:
            #  1. the original subsetter, up to but excluding the insertion pt
            #  2. the insert
            #  3. the remainder of the original subsetter
            new_subsetter = (self._subsetters[k][:subsetter_idx] +
                subsetter_insert + self._subsetters[k][subsetter_idx:])
            self._subsetters[k] = new_subsetter


    # draw recombination rates as needed according to parameter values
    def _draw_recombination_rates(self):
        if (self._r_distr_alpha is not None and self._r_distr_beta is not None):
            recomb_rates = np.clip(np.random.beta(a=self._r_distr_alpha,
                                                  b=self._r_distr_beta,
                                                  size=len(self._positions)),
                                   a_min=0, a_max=0.5)
        elif self._r_distr_alpha is not None:
            recomb_rates = np.ones(len(self._positions)) * self._r_distr_alpha
        else:
            #NOTE: if no alpha and beta values were provided for the recomb-rate
            # beta distribution, set all interlocus recomb rates to a value that
            # results in an average of one recombination per gamete produced
            homog_recomb_rate = 1 / self._L
            recomb_rates = np.ones(len(self._positions)) * homog_recomb_rate
        # set the 0th rate to 0 (because this will ensure that all recombination
        # paths start on homologue 0, such that when the start homologue of a
        # certain gamete equals 1 and its genome array is L-R flipped before
        # subsetting, the effective start homologue of the new, subsetted gamete
        # will actually be 1)
        recomb_rates[0] = 0
        return recomb_rates


    # simulate recombination events
    def _draw_recombination_events(self, use_subsetters, nonneutral_loci,
                                   use_tskit):
        """
        NOTE: Positions and recomb_rates must be provided already sorted!
        """
        # create the breakpoints dict
        recombinations = [r.binomial(1, self._rates) for _ in range(self._n)]
        breakpoints = [self._positions[np.where(
                                        recomb)] for recomb in recombinations]
        for bp in breakpoints:
            if len(bp) > 0:
                assert 0 not in bp, '%s' % str(bp)
        # recast it as an int-keyed dict, so that when I use the recombination
        # events to simulate recombination I don't have to pass around the long
        # arrays, but instead just random keys to index them out on the fly
        breakpoints = dict(zip(range(len(breakpoints)), breakpoints))

        if use_subsetters:
            # get the recombination paths used to make the bitarray subsetters
            if ((use_tskit and len(nonneutral_loci) > 0)
                or not use_tskit):
                recomb_paths = [np.cumsum(
                                    recomb) % 2 for recomb in recombinations]
                # if using tskit, grab just the values at the nonneutral loci
                # (because the subsetters will be used only for subsetting
                # the non-neutral genotypes carried with the individuals)
                if use_tskit:
                    paths = [path[nonneutral_loci] for path in recomb_paths]
                # otherwise, use all values
                else:
                    paths = [path for path in recomb_paths]

                subsetters = [bitarray.bitarray(''.join([self._bitarray_dict[
                            step] for step in path])) for path in paths]
                # recast as integer-indexed dict
            else:
                subsetters = [bitarray.bitarray() for _ in range(
                                                        len(recombinations))]
            subsetters = dict(zip(range(len(subsetters)), subsetters))
            assert len(breakpoints) == len(subsetters)
        else:
            subsetters = None
        return (breakpoints, subsetters)


    # calculate the left and right segment edges for each recombination event
    def _set_seg_info(self):
        seg_info = {}
        ##@##
        for k, recomb in self._breakpoints.items():
            # NOTE: subtracting 0.5 from all recombination breakpoints,
            # to indicate that recombination 'actually' happens halfway between
            # a pair of loci, (i.e. it actually subsets Individuals'
            # genomes in that way) without having the hold all the 0.5's
            # in the Recombinations._events data struct (to save on memory)
            left = np.array([0] + [*recomb - 0.5])
            right = np.array([*recomb - 0.5] + [self._L])
            assert np.all((right-left)>0), ("Some of the following breakpoints"
                        " are invalid: :%s") % str([*zip(left, right)])
            seg_info[k] = (left, right)
        self._seg_info = seg_info


   # take a recombination event key and a parent's node ids,
   # return a zip object containing 
        # 1.) node id for parent's id (corresponding to the id in the
        # tskit.TableCollection.nodes table);
        # 2.) left end of the segment (inclusive, according to tskit)
        # 3.) right end (exclusive)
    def _get_seg_info(self, start_homologue, event_key,
                                    node_ids):
        left, right = self._seg_info[event_key]
        # add small jitter to the numbers, so that recombination breakpoints
        # are essentially not repeated, and this way tskit can correctly track
        # the number of breakpoints and number of trees
        if self._jitter_breakpoints:
            left[1:] = left[1:] + np.random.uniform(0, 0.0001, len(left)-1)
            #left, right = [np.clip(seg_bps + np.random.uniform(0, 0.0001,
            #                                                   len(seg_bps)),
            #                       0, self._L) for seg_bps in [left, right]]

            # NOTE: must make sure that each set of seg_info
            # always starts at 0 on the leftmost segment
            # (and ends at gen_arch.L on the rightmost, but that will already
            #  happen in the above code because we're only ever adding a
            #  positive value to original rightmost value, which is gen_arch.L,
            #  then clipping that value back down to gen_arch.L)
            #left[0] = 0
            right[:-1] = left[1:]

        homologue_nodes = node_ids[[(i + start_homologue) % 2 for i in range(
                                                                   len(left))]]
        seg_info = zip(homologue_nodes, left, right)
        return seg_info


class Trait:
    """
    Representation of a single trait and its parameterization.

    Trait objects are collected in a serial integer-keyed dict,
    one per trait that is parameterized for the corresponding Species,
    which inheres to the Species' GenomicArchitecture object as the attribute
    `GenomicArchitecture.traits`, from which individual Traits can be indexed
    using their index-number keys (e.g. `GenomicArchitecture.traits[<idx>]`).

    Attributes
    ----------

        NOTE: For more detail, see the documentation for the parameters
              that correspond to many of the following attributes.

        alpha:
            A 1d numpy array containing the effect sizes of all loci subtending
            this trait. The value at index i represents the effect size of the
            locus at index i in the 'loci' attribute of this Trait.

        alpha_distr_mu:
            Mean of the normal (Gaussian) distribution
            of effect sizes for this trait's loci.

        alpha_distr_sigma:
            Standard deviation of the normal (Gaussian) distribution
            of effect sizes for this trait's loci.

        gamma:
            Curvature of the fitness function.

        idx:
            Index number of the trait (i.e. its key within the traits dict)

        loci_idxs:
            A 1d numpy array containing the index positions within an
            Individual's genotypes array
            (attribute 'g' of an Individual) that
            correspond to the loci subtending this trait.
            (This is designed such that an Individual's genotype
            for the locus number stored at Trait.loci[i]
            can be accessed by calling Individual.g[Trait.loci_idxs[i], :].)
            This method of tracking allows Individuals' full genotypes to be
            stored in tskit tables (if tskit is being used)
            while they carry their non-neutral genotypes with them
            (a computational optimization).
            Set to None if tskit is not being used.

        loci:
            A numerically sorted 1d numpy array containing
            the locus numbers of all loci subtending this trait.

        lyr_num:
            Index number of the Landscape layer to which this trait corresponds
            (i.e. of the environmental variable that serves as the selective
            force on this trait)

        max_alpha_mag:
            Maximum absolute value that can be drawn for a locus’
            effect size (i.e. alpha). Effect sizes are clipped
            to the closed interval [-max_alpha_mag, max_alpha_mag]
            when the model is parameterized.

        mu:
            Mutation rate for this trait (expressed in mutations per locus per
            time step)

        name:
            The string name of the trait

        n_loci:
            Number of loci subtending this trait

        phi:
            Polygenic selection coefficient (i.e the selection coefficient
            acting on the phenotypes, rather than the genotypes, of this Trait)

        univ_adv:
            Whether or not this trait will be treated as universally
            advantageous. When False (the default behavior), individuals will
            be fitter when their phenotypic values for this trait more closely
            align with the local environmental value of the landscape layer
            corresponding to this trait (i.e. selection will be spatially
            varying). When False, phenotypes closer to 1 will be more
            advantageous globally on the Landscape.
    """
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
        self.loci_idxs = np.int64([])
        self.alpha = np.array([])

    def _get_phi(self, spp):
        if type(self.phi) in (float, int):
            phi = np.array([self.phi]*len(spp))
        else:
            phi = self.phi[spp._cells[:, 1], spp._cells[:, 0]]
        return(phi)

    def _set_loci(self, loci):
        # set the loci
        self.loci = np.hstack((self.loci, np.array([*loci])))
        self.loci.sort()
        # set the number of loci
        self.n_loci = self.loci.size

    def _set_loci_idxs(self, nonneut_loci, use_tskit):
        # if using tskit then only non-neutral loci will be stored in
        # Individuals' genotype arrays,
        # so set an index of the loci's indices in the nonneut_loci object,
        # to use when subsetting for the trait's genotypes
        if use_tskit:
            self.loci_idxs = np.array([np.where(
                                nonneut_loci == n)[0][0] for n in self.loci])
        else:
            self.loci_idxs = None

    def _add_locus(self, locus, alpha, idx):
        # get insertion index for the locus 
        insert_pt = bisect.bisect_left(self.loci, locus)
        # insert in the loci
        self.loci = np.hstack((self.loci[:insert_pt],
                               locus,
                               self.loci[insert_pt:]))

        # insert the effect size for the locus
        self.alpha= np.hstack((self.alpha[:insert_pt],
                               alpha,
                               self.alpha[insert_pt:]))

        # insert in the loci_idxs
        # NOTE: add 1 to all superseding locus indexes, to account for
        # locus newly inserted into genotype arrays
        self.loci_idxs = np.hstack((self.loci_idxs[:insert_pt],
                                  idx,
                                  self.loci_idxs[insert_pt:] + 1))

        # increment the number of loci
        self.n_loci += 1


class GenomicArchitecture:
    """
    Representation of the genomic architecture of a Species.

    Inheres to its Species as the Species.gen_arch attribute.

    Attributes
    ----------

        NOTE: For more detail, see the documentation for the parameters
              that correspond to many of the following attributes.

        delet_alpha_distr_scale:
            The scale parameter of the Gamma distribution from which
            the effect sizes (i.e. alphas) of deleterious loci are drawn

        delet_alpha_distr_shape:
            The shape parameter of the Gamma distribution from which
            the effect sizes (i.e. alphas) of deleterious loci are drawn

        delet_loci:
            A 1d numpy array that tracks all deleterious loci

        delet_loci_idxs:
            A 1d numpy array containing the index positions within an
            Individual's genotypes array (attribute 'g' of an Individual)
            that correspond to the deleterious loci.
            (This is designed such that an Individual's genotype for
            the deleterious locus number stored at
            GenomicArchitecture.delet_loci[i] can be accessed by calling
            Individual.g[GenomicArchitecture.delet_loci_idxs[i], :].)
            This method of tracking allows Individuals' full genotypes to be
            stored in the tskit tables (if tskit is being used)
            while they carry their non-neutral genotypes with them
            (a computational optimization).
            Set to None if tskit is not being used.

        delet_loci_s:
            A 1d numpy array containing the selection strengths
            of each of the deleterious loci. (Selection strengths are
            listed in an order corresponding to the order of
            GenomicArchitecture.delet_loci, such that the selection strength of
            the deleterious locus referenced at
            GenomicArchitecture.delet_loci[i] is stored at
            GenomicArchitetcture.delet_loci_s[i].)

        dom:
            Dominance values for all loci (stored as an L-length,
            1d numpy array in which 0 indicates a codominant locus
            and 1 indicates a locus for which the '1' allele is dominant)

        L:
            Length of the simulated genome

        mu_delet:
            Genome-wide deleterious mutation rate (expressed in mutations
            per locus per time step). Deleterious mutations do not influence
            trait phenotypes, but instead are univerally
            deleterious (to an extent determined by their effect sizes).

        mu_neut:
            Genome-wide neutral mutation rate (expressed in mutations
            per locus per time step).

        neut_loci:
            A 1d numpy array that tracks all loci that do not influence the
            phenotypes of any traits.

        nonneut_loci:
            A 1d numpy array that tracks all loci that either influence the
            phenotype of at least one trait or are deleterious.

        p:
            Starting allele frequencies for all loci (stored as an L-length,
            1d numpy array)

        pleiotropy:
            A bool indicating whether or not to allow pleiotropy (i.e. whether
            or not to allow the same locus to subtend multiple traits)

        recombinations:
            The Recombinations object, which contains the parameters and data
            necessary for simulation of recombination. Not intended for public
            use at this time.

        sex:
            A bool indicating whether or not Individuals of this Species
            will be assigned a sex (0 for female, 1 for male)

        traits:
            A dict containing all of a Species' Trait objects, keyed by serial
            integers.

        x:
            Ploidy
            NOTE: only diploidy is currently implemented, but this attribute
                  serves as a placeholder for use in possible
                  future generalization to arbitrary ploidy.)
    """
    def __init__(self, dom, g_params, land, recomb_rates=None,
                 recomb_positions=None):
        # ploidy (NOTE: for now will be 2 by default; later could consider
        # enabling polyploidy)
        self.x = 2
        # total length (i.e. number of markers)
        self.L = g_params.L
        # placeholder for the starting allele freqs for all loci
        self.p = None
        # True/False regarding whether to allow a locus to affect the
        # phenotype of more than one trait; defaults to False
        self.pleiotropy = g_params.pleiotropy
        # array of dominance values for all loci
        self.dom = dom
        # set the _use_dom attribute based on whether any loci have a 1 for
        # their dominance value
        self._use_dom = np.any(self.dom)
        # whether or not to use sexes for this species
        self.sex = g_params.sex

        # whether or not to use tskit to track this species' spatial pedigree
        self.use_tskit = g_params.use_tskit
        assert self.use_tskit in [True, False], (("'use_tskit' must be either "
                                                  "True or False for each "
                                                  "Species."))
        # set interval for simplifying the tskit TableCollection
        self.tskit_simp_interval = g_params.tskit_simp_interval
        if self.use_tskit:
            assert isinstance(self.tskit_simp_interval,
                              int), (("If using tskit then the "
                                      "'tskit_simp_interval' parameter must be "
                                      " an integer."))

        # genome-wide neutral mutation rate
        self.mu_neut = g_params.mu_neut
        # array to keep track of all loci that don't influence the
        # phenotypes of any trait; defaults to all loci, then will be updated
        self.neut_loci = np.array(range(self.L))
        # array to keep track of all loci that influence the phenotype of at
        # least one trait; after burn-in, will be updated as needed
        self.nonneut_loci = np.array([])

        # genome-wide deleterious mutation rate
        self.mu_delet = g_params.mu_delet
        self.delet_alpha_distr_shape = g_params.delet_alpha_distr_shape
        self.delet_alpha_distr_scale = g_params.delet_alpha_distr_scale


        # arrays to track deleterious loci, their genotype indices, and their
        # strengths of selection
        self.delet_loci = np.int64([])
        self.delet_loci_idxs = np.int64([])
        self.delet_loci_s = np.array([])

        # add a dict of Trait objects, if necessary
        self.traits = None
        if 'traits' in [*g_params]:
            self.traits = _make_traits(g_params.traits, land)

        # set self._mu_tot, the total per-site, per-generation mutation rate
        mus = [mu for mu in (self.mu_neut, self.mu_delet) if mu is not None]
        if self.traits is not None:
            mus = mus + [trt.mu for trt in self.traits.values()]
        self._mu_tot = sum(mus)
        # save the nonneutral mutation rate
        self._mu_nonneut = self._mu_tot - self.mu_neut

        # set a placeholder for the species' mutable loci
        # (to be filled after burn-in)
        self._mutables = None

        # attribute that will be replaced with the estimated total number of
        # mutations per iteration, once it's estimated at the end of the first
        # burn-in
        self._mut_fns = self._make_mut_fns_dict()
        # set ._planned_muts to None, for now (this is not yet implemented,
        # but thinking about it
        self._planned_muts = None

        # The recombination-paths object will be assigned here; used to
        # speed up large quantities of binomial draws needed for recombination
        self.recombinations = Recombinations(self.L, recomb_positions,
                                             g_params.n_recomb_sims,
                                             g_params.r_distr_alpha,
                                             g_params.r_distr_beta,
                                             recomb_rates,
                                             g_params.jitter_breakpoints)

    # method to make a _mut_fns dict, containing a function
    # for each type of mutation for this species
    def _make_mut_fns_dict(self):
        mut_fns = {}
        if self.mu_neut > 0:
            def neut_fn(spp, offspring):
                return(mutation._do_neutral_mutation(spp, offspring))
            mut_fns.update({'neut': neut_fn})
        if self.mu_delet > 0:
            def delet_fn(spp, offspring):
                return(mutation._do_deleterious_mutation(spp, offspring))
            mut_fns.update({'delet': delet_fn})
        if self.traits is not None:
            for trait_num in self.traits:
                if self.traits[trait_num].mu > 0:
                    def trait_fn(spp, offspring, trait_num=trait_num):
                        return(mutation._do_trait_mutation(spp, offspring,
                                                           [trait_num]))
                    mut_fns.update({'t%i' % trait_num: trait_fn})
        return mut_fns

    # method to draw mutation types for any number of mutations chosen
    # to occur in a given timestep
    def _draw_mut_types(self, num):
        type_dict = {'neut': self.mu_neut,
                     'delet': self.mu_delet}
        if self.traits is not None:
            trait_dict = {'t%i' % (k): v.mu for k, v in self.traits.items()}
            type_dict.update(trait_dict)
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
        if not mutational:
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
            loci = set(r.choice(self.neut_loci, size=n, replace=False))
        elif self.pleiotropy:
            loci = set(r.choice(range(self.L), size=n, replace=False))
        # update the trait's loci
        self.traits[trait_num]._set_loci(loci)
        # add these loci to self.non-neutral and remove from
        # self.neut_loci, to keep track of all loci underlying any
        # traits (for purpose of avoiding pleiotropy)
        self.nonneut_loci = np.array(sorted([*self.nonneut_loci] + [*loci]))
        self.neut_loci = np.array(sorted([*set(self.neut_loci).difference(
                                                    set(self.nonneut_loci))]))
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

    # add a nonneutral locus to the genomic architecture
    # NOTE: either trait_nums or delet_s must be non-None,
    # and trait_nums must be iterable (if not None)
    def _add_nonneut_locus(self, locus, trait_nums=None, delet_s=None):
        # remove from the neut_loci array
        self.neut_loci = np.delete(self.neut_loci,
                                   np.where(self.neut_loci == locus))
        # add the locus to the nonneut_loci array
        idx = bisect.bisect_left(self.nonneut_loci, locus)
        self.nonneut_loci = np.hstack((self.nonneut_loci[:idx],
                                       locus,
                                       self.nonneut_loci[idx:]))

        # add the locus to either trait loci, if necessary
        if trait_nums is not None and delet_s is None:
            for n in trait_nums:
                alpha = self._draw_trait_alpha(n)[0]
                self.traits[n]._add_locus(locus, alpha, idx)
        # or else add to the deleterious loci
        elif delet_s is not None and trait_nums is None:
            # TODO: DECIDE IF NEED TO USE self.loci INSTEAD OF self.delet_loci
            #       HERE IF use_tskit==False
            del_idx = bisect.bisect_left(self.delet_loci, locus)
            # add the locus, its genome-index (for subsetting individuals'
            # genomes when calculating fitness), and its strength of selection
            # to the deleterious locus trackers
            self.delet_loci = np.hstack((self.delet_loci[:del_idx],
                                         locus,
                                         self.delet_loci[del_idx:]))
            if self.use_tskit:
                self.delet_loci_idxs = np.hstack(
                    (self.delet_loci_idxs[:del_idx], idx,
                     self.delet_loci_idxs[del_idx:]))
            else:
                self.delet_loci_idxs = None
            self.delet_loci_s = np.hstack((self.delet_loci_s[:del_idx],
                                         delet_s,
                                         self.delet_loci_s[del_idx:]))
        return idx


    # method for plotting all allele frequencies for the species
    def _plot_allele_frequencies(self, spp, color='red'):
        if spp.gen_arch.use_tskit:
            speciome = spp._get_genotypes(biallelic=False)
        else:
            speciome = np.mean(np.stack([ind.g for ind in spp.values()]), axis=2)
        freqs = np.mean(speciome, axis=0)
        assert len(freqs) == spp.gen_arch.L
        plt.plot(range(self.L), self.p, ':k', label='start freq.', alpha=0.5)
        plt.plot(range(self.L), freqs, '-', color=color, label='curr. freq.')
        plt.xlabel('locus')
        plt.ylabel('frequency')
        plt.legend()
        plt.show()

    # method for pickling a genomic architecture
    def _write_pickle(self, filename):
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


# build the genomic architecture
def _make_genomic_architecture(spp_params, land):
    # get the genome parameters
    g_params = spp_params.gen_arch
    # get the custom genomic-architecture file, if provided
    gen_arch_file = None
    if 'gen_arch_file' in g_params.keys():
        if g_params.gen_arch_file is not None:
            gen_arch_file = pd.read_csv(g_params.gen_arch_file)
            # check that the genome length matches the params file
            assert len(gen_arch_file) == g_params.L, ("The length of the "
                        "custom genomic architecture file must match the "
                        "genome length provided to the parameter 'L' in "
                        "Geonomics parameters file.")
            # check that each trait is assigned an alpha value
            for trt_val, alpha_val in zip(gen_arch_file['trait'],
                                          gen_arch_file['alpha']):
                if pd.notnull(trt_val):
                    assert pd.notnull(alpha_val), ("All trait-associated loci "
                                                   "specified in a custom "
                                                   "genomic architecture file "
                                                   "must be assigned non-null "
                                                   "alpha values. Please check "
                                                   "your file, then rerun.")

            # ensure that trait and alpha columns are strings to start (because
            # in many cases alpha would likely read in as a float, if each row
            # has at most one alpha value because there's no pleiotropy, but in
            # case of pleiotropy I need to be able to use the str.split method
            # on each row's values further down)
            gen_arch_file['trait'] = [str(v) for v in gen_arch_file['trait']]
            gen_arch_file['alpha'] = [str(v) for v in gen_arch_file['alpha']]
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

    # get the sex parameter and add it as an item in g_params
    g_params['sex'] = spp_params.mating.sex

    # get the custom recomb_rates and recomb_positions from the
    # custom gen-arch file, if provided
    # NOTE: each recombination rate
    #       stipulated in gen_arch file
    #       is construed as the recombination rate
    #       at the midway point between that locus
    #       and the previous one
    #       (NOTE: this is of course non-sensical for the 0th
    #              locus, which would thus have its rate expressed as
    #              the rate of recomb between the -1th and the 0th,
    #              but that's fine because the first rate is later coerced to 0
    #              anyhow, then the starting homologue of the recomb path,
    #              i.e. the sister chromatid placed into the gamete,
    #              is later chosen using a ~Bern(0.5) draw at the time
    #              mating occurs; see mating.py)
    recomb_rates = None
    recomb_positions = None
    if gen_arch_file is not None:
        recomb_rates = gen_arch_file['r'].values
        recomb_positions = gen_arch_file['locus'].values

    # set locus-wise dominance values for the 1-alleles, using the 'dom' value
    # in the gen_arch params, unless a gen_arch_file was provided
    if gen_arch_file is None:
        # create an L-length array of boolean integers (where 0 = codominance,
        # 1 = dominance)
        dom = np.array([int(g_params.dom)] * g_params.L)
    else:
        # get the 'dom' column from the gen_arch_file
        dom = gen_arch_file['dom'].values

    # now make the gen_arch object
    gen_arch = GenomicArchitecture(dom, g_params, land,
                                   recomb_rates, recomb_positions)

    # set the loci and effect sizes for each trait, using the custom gen-arch
    # file, if provided, and if the species has traits
    if gen_arch_file is not None and gen_arch.traits is not None:
        # convert the trait names in the 'trait' column of the file into
        # lists of their trait numbers (i.e. their keys in the gen_arch
        # traits dict)
        trt_names_nums = {
            trt.name: num for num, trt in gen_arch.traits.items()}
        gen_arch_file['trait'] = [[trt_names_nums[
            val] for val in [x.strip() for x in row.split(
            ',')] if val in trt_names_nums] for row in gen_arch_file['trait']]
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
        # check that we got the same length of loci and effect sizes for
        # a trait, for all traits
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

    # now that all nonneutral loci have been selected, 
    # set each trait's loci_idxs (which maps a trait's numeric loci to
    # their index positions in individuals' genomes)
    # NOTE: will just be set to None if not using tskit
    if gen_arch.traits is not None:
        for trt in gen_arch.traits.values():
            trt._set_loci_idxs(gen_arch.nonneut_loci, gen_arch.use_tskit)

    assert_msg = ("The union of the gen_arch.neut_loci and "
                  "gen_arch.nonneut_loci sets does not contain all loci "
                  "indicated by gen_arch.L")
    assert len(set(range(gen_arch.L)).difference(
        set(gen_arch.neut_loci).union(
        set(gen_arch.nonneut_loci)))) == 0, assert_msg

    # draw locus-wise 1-allele frequencies, unless provided in
    # custom gen-arch file
    if gen_arch_file is None:
        # set according to a non-None value
        if g_params.start_p_fixed is not None:
            if isinstance(g_params.start_p_fixed, bool):
                # if True, make all 0.5 (default starting allele freq)
                if g_params.start_p_fixed:
                    gen_arch.p = np.array([0.5]*g_params.L)
                # if False, draw randomly
                else:
                    gen_arch.p = _draw_allele_freqs(g_params.L)
            # if not a bool then should be a float, so assert,
            # then set to that val
            else:
                assert 0 <= g_params.start_p_fixed <= 1, ("If a starting allele"
                                                        " frequency value is "
                                                        "provided then it must "
                                                        "be between 0 and 1.")
                gen_arch.p = np.array([g_params.start_p_fixed]*g_params.L)
        # if the value is None, then draw randomly
        else:
            gen_arch.p = _draw_allele_freqs(g_params.L)
        # set the neutral loci to 0, if need be
        if g_params.start_neut_zero and len(gen_arch.neut_loci) > 0:
            gen_arch.p[gen_arch.neut_loci] = 0
    # use starting allele freq values from the gen_arch_file instead, if given
    else:
        gen_arch.p = gen_arch_file['p'].values

    # set the recombination events
    use_subsetters = ((gen_arch.use_tskit and (len(gen_arch.nonneut_loci) > 0 or
                                               gen_arch._mu_nonneut > 0))
                      or not gen_arch.use_tskit)
    gen_arch.recombinations._set_events(use_subsetters, gen_arch.nonneut_loci,
                                        gen_arch.use_tskit)

    return gen_arch


# reset genomes after burn-in
def _check_mutation_rates(gen_arch, est_tot_muts, burn_T, T):
    # if there do not appear to be enough space in the simulated genome for
    # the expected number of mutations then warn the user of that
    if est_tot_muts > 0.75 * len(
        [loc for loc in range(
                    gen_arch.L) if loc not in gen_arch.nonneut_loci]):
        raise MutationRateError(("This species has been "
                                 "parameterized with too few neutral "
                                 "loci to accommodate the expected "
                                 "number of mutations. (Geonomics only "
                                 "uses an infinite sites model.) "
                                 "Please tweak some combination of "
                                 "the genome length, model run time, "
                                 "or mutation rates."))


    # skip this step and force the neutral mutation rate to 0, if there are no
    # neutral loci in the genome as it was configured
    if len(gen_arch.neut_loci) == 0 and gen_arch._mu_tot > 0:
        warn_msg = ("This species has been parameterized "
                    "with non-zero mutation rates but "
                    "without any neutral loci, leaving no target "
                    "for mutations. Please tweak the genome "
                    "length and/or mutation rates.")
        warnings.warn(warn_msg)
        gen_arch.mu_neut = 0
        gen_arch.mu_delet = 0
        for trt in gen_arch.traits.values():
            trt.mu = 0
    # or just outright skip this step if the parameterization included only
    # mutation rates of 0
    elif gen_arch._mu_tot == 0:
        pass
    # otherwise, set the mutable loci
    else:
        mutables = [*set(range(gen_arch.L)).difference(
                                            set(gen_arch.nonneut_loci))]
        r.shuffle(mutables)
        gen_arch._mutables = [*mutables]
    return


def _make_starting_mutations(spp, tables):
    '''
    function to generate mutations, after burn-in,
    and to assign them to a species' TableCollection's current nodes,
    to produce the starting 1-allele frequencies parameterized for the species
    '''
    # get the starting frequencies for each site
    start_freqs = spp.gen_arch.p

    # get a set of all homologues, as tuples of (individual id, homologue idx)
    homologues = [*zip(np.repeat([*spp], 2),
                         [*range(spp.gen_arch.x)] * len(spp))]

    # make mutations for each site
    for site, freq in enumerate(start_freqs):
        # generate the number of mutations for this locus
        n_mutations = int(round(2 * len(spp) * freq, 0))
        # make sure we don't mutate either all or none of the population's
        # homologues, unless called for
        if n_mutations == len(spp) * 2 and freq < 1:
            n_mutations -= 1
        if n_mutations == 0 and freq > 0:
            n_mutations = 1
        #randomly choose and mutate n_mutations homologues from the population 
        np.random.shuffle(homologues)
        homologues_to_mutate = homologues[:n_mutations]
        for ind, homol in homologues_to_mutate:
            # if using tskit, add muts to individis' genotype arrays as needed,
            # and also to mutations table
            if spp.gen_arch.use_tskit:
                # create a mutation in the individual's genome, if this is a
                # non-neutral locus
                if site in spp.gen_arch.nonneut_loci:
                    spp[ind].g[np.where(spp.gen_arch.nonneut_loci == site),
                               homol] = 1
                # get the homologue's nodes-table id,
                # then add a row to the mutations table
                node_id = spp[ind]._nodes_tab_ids[homol]
                tables.mutations.add_row(site, node=node_id, parent=-1,
                                         derived_state='1')

            # if not using tskit just add all muts to individs' genotypes arrays
            else:
                spp[ind].g[site, homol] = 1

    # and then reset the individuals' phenotypes, if needed
    if spp.gen_arch.traits is not None:
        [ind._set_z(spp.gen_arch) for ind in spp.values()]

    return


# method for loading a pickled genomic architecture
def read_pickled_genomic_architecture(filename):
    import cPickle
    with open(filename, 'rb') as f:
        gen_arch = cPickle.load(f)
    return gen_arch


################################################################
# FUNCTIONS FOR WORKING WITH msprime TO SIMULATE GNX INDIVIDUALS
################################################################

def _format_msprime_rate(rate, L):
    """
    Either return a simple float, if the given rate
    is homogeneous, or else return an msprime.RateMap object.
    Used for msprime recombination and mutation rates.
    """
    assert (type(rate) in [float, list, np.ndarray] or rate == 0 or rate == 1)
    if (isinstance(rate, float) or rate == 0 or rate == 1):
        return rate
    elif len(np.unique(rate)) == 1:
        return rate[0]
    else:
        assert len(rate) == L
        # NOTE: position array has to be one longer than the sequence length,
        #       so that rightmost position brackets a stretch of genomic
        #       intervals equal to the sequence_length
        positions = np.linspace(0, L, L+1)
        ratemap = msprime.RateMap(position=positions, rate=rate)
        return ratemap


def _sim_msprime_individuals(n,
                             L,
                             recomb_rate,
                             mut_rate,
                             demography=None,
                             population_size=None,
                             ancestry_model=None,
                             random_seed=None,
                             use_tskit=True,
                            ):
    '''
    Use an msprime simulation to generate new Geonomics Individuals.
    (docs in the main.py wrapper)
    '''
    # handle the recombination and mutation rates
    if (not isinstance(recomb_rate, float) and
        not (recomb_rate == 0)):
        # NOTE: tack a 0.5 onto the righthand side of the recomb_rate vector,
        #       to satisfy the difference between msprime (for which L
        #       expresses length in intervals, each of which needs a recomb
        #       rate specified, hence the recomb_rate vector needing to have a
        #       number of rates equal to L) and gnx (for which L expressed
        #       length in loci, between which recomb rates are needed, hence
        #       needing L-1 recomb rates);
        #       TODO: DOUBLE-CHECK THE FOLLOWING LOGIC
        #       this 0.5 value will not affect gnx data because it effectively
        #       sits between sites L-1 and L, but a simulated gnx genome of
        #       length L only includes sites [0, L-1]
        recomb_rate = [*recomb_rate] + [0.5]
    recomb_rate = _format_msprime_rate(recomb_rate, L)
    mut_rate = _format_msprime_rate(mut_rate, L)

    # simulate ancestry and mutations
    # TODO: look back at how gnx uses genome positions, because 
    #       it looks like we need to provide discrete_genome=false
    #       to sim_ancestry and true to sim_mutations;
    #       or perhaps instead just fine to use false for both
    #       and then subtract 0.5 from all non-zero 'left' values
    #       from all non-L 'right' values when handles the edges table?
    ts = msprime.sim_ancestry(samples=n,
                              sequence_length=L,
                              demography=demography,
                              population_size=population_size,
                              recombination_rate=recomb_rate,
                              ploidy=2,
                              model=ancestry_model,
                              random_seed=random_seed,
                              discrete_genome=True,
                             )
    # NOTE: mutations use a simple binary 0|1 model, with 0 always the
    #       ancestral state, and with mutations being state-independent
    #       (i.e., always placing a 1, even if this results in a silent
    #       mutation of 1→1);
    #       this model is the best aligned to gnx data
    mut_model = msprime.BinaryMutationModel(state_independent=True)
    mts = msprime.sim_mutations(tree_sequence=ts,
                                rate=mut_rate,
                                model=mut_model,
                               )
    # get the mutated tskit.TableCollection object
    mtc = mts.dump_tables()

    # double-check the number of Individuals
    assert np.sum(mtc.nodes.asdict()['individual']>-1) == 2*n

    # create the Individuals, assigning placeholder values for the attributes
    # (to be later updated by downstream code)
    individs = [Individual(idx=0, x=0, y=0) for i in range(n)]
    # get the information linking Individuals to rows in TableCollection tables
    nodes_tab_individs_col = mtc.nodes.asdict()['individual']
    for i, ind in enumerate(individs):
        # get the nodes table rows corresponding to this Individual
        nodes_tab_ids = np.where(nodes_tab_individs_col==i)[0]
        assert len(nodes_tab_ids) == 2
        # set hidden attributes tracking the individs' rows in the 'nodes' table
        # and id in the 'individuals' table
        ind._set_nodes_tab_ids(nodes_tab_ids)
        ind._individuals_tab_id = i
        # if the gnx Model these Individuals are intended for is not using tskit
        # then walk the variants generator and compile the individ's
        # genotype array, then store it in the attribute 'g'
        if not use_tskit:
            var_gen = mts.variants(samples=nodes_tab_ids)
            gt = np.stack([next(var_gen).genotypes for _ in range(L)])
            ind.g = gt

    return individs, mtc


###############################################################################
# FUNCTIONS FOR WORKING WITH THE SPATIAL PEDGREE SAVED IN tskit DATA STRUCTURES
###############################################################################

def _prep_tskit_tabcoll_for_gnx_spp_union(tc,
                                          recip_spp,
                                          coords,
                                          source_software,
                                         ):
    '''
    Reformat a tskit.TableCollection, produced by either a gnx simulation
    or an msprime simulation, so that it conforms to Geonomics' expected format
    and has the other data it needs to be fed into the .union() method
    of the given Geonomics Species' current TableCollection.

    The code and comments walk through each table, column by column,
    and specify if and how the contents need to be reformatted or updated
    in order to match Geonomics conventions.

    Returns the list of gnx Individual idx numbers that were set in the
    TableCollection, for downstream use and/or validation.
    '''
    # check source_software is valid
    assert source_software in ['gnx', 'msprime'], ('This function can only '
                                                   'prep TableCollections '
                                                   'that originate from '
                                                   'gnx Models or msprime '
                                                   'simulations.')

    # use the coords array to get the number of new individuals being unioned
    # into the recipient population
    assert (isinstance(coords, np.ndarray) and
            len(coords.shape)==2 and
            coords.shape[1] == 2), ("The coords param must be a 2d "
                                    "numpy.ndarray with the x,y coords in "
                                    "the 2 columns.")
    n_source_individs = coords.shape[0]
    if source_software == 'msprime':
        assert tc.individuals.num_rows == n_source_individs

    # reformat nodes table content:
    #-----------------------------

    if source_software == 'msprime':
        # gnx uses the time column as follows:
            # - smallest values are the most recent individuals
            # - nodes with val 1 (i.e. -1 * t=-1 before main phase begins) are
            #   nodes that originated at the end of the burn-in phase
            # - nodes with val 0, -1, -2, ... (i.e., -1 * t =0, 1, 2,...)
            #   originated during those time steps
        # msprime times are expressed as t before present (i.e. samples are 0s).
        # thus, convert msprime to gnx convention by subtracting from the whole
        # column the value of spp.t at the time of their introduction
        # (i.e., the time of their simulation and sampling)
        time = tc.nodes.time - recip_spp.t
        assert np.all([t>=(-1*recip_spp.t) for t in time]), ("Values "
                            "attempting to be set on 'time' column of nodes "
                            "table are invalid (i.e., are more negative than "
                            "-1 * recipient Species' current model time, and "
                            "thus would refer to future time steps not yet "
                            "simulated.")
        tc.nodes.time = time

        # msprime uses the flags column to indicate extant nodes
        # (i.e., docs say "NODE_IS_SAMPLE = 1"); gnx uses flags to distinguish
        # between nodes in 'real' gnx Individuals and msprime nodes;
        # thus, we actually just want to keep this column as is, because the
        # samples are only nodes that will be assigned to new gnx Individuals,
        # and thus are the only ones we would otherwise coerce to 1s anyhow!

        # nothing to do the population column; TableCollection.union()
        # will be asked to use the next available int in the recipient Species'
        # populations table to automatically record the
        # population these nodes came from

        # the individual column needn't be touched; TableCollection.union()
        # will automatically increment this col by using the next individual id
        # available in the recipient Species' individuals table;
        # we will just need to later track what the next id was and how many
        # individuals were added, so that TableCollection.individuals ids can
        # be matched to gnx Individuals via Individuals._individuals_tab_id attr

        # metadata, metadata_offset, and metadata_schema columns are currently
        # unused by both msprime and gnx

    elif source_software == 'gnx':
        # no need to alter any of this data if it's already from a gnx model
        pass



    # reformat individuals table content:
    #-----------------------------------

    if source_software == 'msprime':
        # all individuals in the individuals table are extant and will be turned
        # into gnx Individual objects, and gnx uses a 1 in the 'flags' column to
        # indicated 'real' Individuals (i.e., individuals that were actual
        # gnx.Individuals at some point), so set the entire 'flags' column to 1s
        flags = tc.individuals.flags*0+1
    elif source_software == 'gnx':
        # flags should already be correct for a gnx model,
        flags = tc.individuals.flags

    # gnx stores individs' starting x,y coordinates in the 'location' column
    # (NOTE: also the phenotype and fitness values, but those are not relevant
    # in the neutral-loci-only msprime simulations permitted by gnx);
    # TODO: think about if/how sweeps would be allowed and
    #       how they would be represented/matched up to genarch
    # location_offset tells the TableCollection to jump 2 values (x, y) per row
    # NOTE: 'ragged columns' actually consist of the column's contents in a
    #       flattened array format and a column of n+1 index values that can be
    #       offset and used to index the start and end positions of each row's
    #       contents within that array
    if source_software == 'msprime':
        location_offset = [int(n) for n in np.linspace(0,
                                                   2*tc.individuals.num_rows,
                                                   tc.individuals.num_rows+1)]
    elif source_software == 'gnx':
        # in a TableCollection from a gnx Model, locations will need to be
        # blanked out for all Individuals, current or past (so that fns that
        # later depend on the stored locations information don't accidentally
        # use data conveying locations on a different landscape from before the
        # introduction of these source Individuals to their new population),
        # and then will need to be reset to the new coords values for only the
        # bottom n rows; thus, need to tack leading zeros onto the
        # location_offset values just generated, making the final vector
        # num_rows+1 long
        location_offset = [int(n) for n in np.linspace(0,
                                                   2*n_source_individs,
                                                   n_source_individs + 1)]
        location_offset = [0] * (tc.individuals.num_rows -
                                 n_source_individs) + location_offset
        assert len(location_offset) == (tc.individuals.num_rows + 1), ("The "
                                "location_offset vector must be one longer "
                                "than the number of rows in the individuals "
                                "table on which it will be set.")

    # parents column is unused by geonomics
    # (wasn't implemented until 21/01/21, after gnx was developed;
    # tskit-dev commit 8ebfd8b2900217dab2b82d6f530c830238037fc8).
    # useful but not necessary to implement;
    # this column is also empty in msprime

    # gnx uses the metadata column to store each individual's
    # gnx Individual idx, for later matching to rows in the individuals table
    # when the table is updated after tskit's simplify algo filters individuals;
    # determine the next n indx values using the Species.max_ind_idx attribute
    gnx_individ_idxs = [*range(recip_spp.max_ind_idx+1,
                       recip_spp.max_ind_idx + n_source_individs + 1)]
    metadata=[idx.to_bytes(length=4,
                           byteorder='little') for idx in gnx_individ_idxs]
    # for msprime, make location_offset vals uniform across new Individs' rows
    if source_software == 'msprime':
        metadata_offset = [int(n) for n in np.linspace(0,
                                                   4*tc.individuals.num_rows,
                                                   tc.individuals.num_rows+1)]
    elif source_software == 'gnx':
        # same as for location_offset above, need to tack leading zeros onto
        # the metadata vector if setting columns in a gnx-generated individuals
        # table
        metadata_offset = [int(n) for n in np.linspace(0,
                                                   4*n_source_individs,
                                                   n_source_individs+1)]
        metadata_offset = [0] * (tc.individuals.num_rows -
                                 n_source_individs) + metadata_offset
        assert len(metadata_offset) == (tc.individuals.num_rows + 1), ("The "
                                "metadata_offset vector must be one longer "
                                "than the number of rows in the individuals "
                                "table on which it will be set.")

    # NOTE the set_columns() method fails because the metadata column, although
    #      supposed to be a binary data type according to the docs, is actually
    #      coded as np.int8 and so fails to convert the 4-byte binary data we
    #      are trying to provide; never ran into this problem in tha past
    #      because the .add_row() method allows the binary data for some reason
    #      and that's the only method I've used before (just adding 1
    #      individual at a time as they're born);
    #      however, somewhat conveniently, we were already planning to change
    #      every column already containing data in this table anyhow, so for
    #      now we are just:
    #           1.) double-checking the length of the original table
    #           2.) clearing the table
    #           3.) manually looping and adding new rows
    #           4.) confirming identical resulting length
    if source_software == 'msprime':
        assert (len(flags) ==
            coords.shape[0] ==
            (len(location_offset)-1) ==
            len(metadata) ==
            (len(metadata_offset)-1)), ("The 'flags', 'location', and "
                                        "'metadata' columns for the "
                                        "individuals table must all "
                                        "be the same length as the coords "
                                        "array's 0th dimension and the "
                                        "'location_offset' and "
                                        "'metadata_offset' columns must be "
                                        "one longer.")
    elif source_software == 'gnx':
        assert (len(flags) ==
                (len(location_offset)-1) ==
                (len(metadata_offset)-1)), ("The 'flags' column must be "
                            "one shorter than the 'location_offset' and "
                            "'metadata_offset' columns.")
        assert len(metadata) == coords.shape[0] == n_source_individs, ("The "
                            "'metadata' column must have one value "
                            "per coordinate pair provided.")
    original_inds_tab_len = tc.individuals.num_rows
    if source_software == 'msprime':
        assert original_inds_tab_len == coords.shape[0], ("Length "
                                             "of individuals table "
                                             "generated by msprime "
                                             "differs from length of coords "
                                             "attempting to be added to it.")
    tc.individuals.clear()
    for i, flag in enumerate(flags):
        # for gnx, add original flags and otherwise empty rows for all but the
        # final rows corresponding to the Individuals being unioned into a
        # recipient Species
        if source_software == 'gnx':
            if i < (len(flags) - n_source_individs):
                row_coords = None
                row_metadata = None
            else:
                i_adjusted = i - (len(flags) - n_source_individs)
                row_coords = coords[i_adjusted, :]
                row_metadata = metadata[i_adjusted]
        elif source_software == 'msprime':
            row_coords = coords[i, :]
            row_metadata = metadata[i]
        tc.individuals.add_row(flags=flag,
                               location=row_coords,
                               parents=None,
                               metadata=row_metadata,
                              )
    assert tc.individuals.num_rows == original_inds_tab_len, ("Reconstructed "
                                            "individuals table is not the same "
                                            "length as the original table "
                                            "generated by msprime.")


    # reformat edges table content:
    #-----------------------------

    # geonomics models all edges as being delineated by positions belonging
    # to the set [0, 0.5, 1.5, ..., L-2.5, L-1.5, L];
    # this is more complicated than the saner msprime
    # convention of using the set [0, 1, 2, ..., L-2, L-1, L],
    # but this was a decision I made based on modeling recombination
    # breakpoints halfway between loci and it is equally functional,
    # so for msprime-derived tables we need to adjust the positions
    if source_software == 'msprime':
        L = recip_spp.gen_arch.L
        left = tc.edges.left
        right = tc.edges.right
        left = np.clip(left - 0.5, a_min=0, a_max=L)
        right = np.clip([r-0.5 if r != L else r for r in right], a_min=0, a_max=L)
        tc.edges.left = left
        tc.edges.right = right

    # parent and child columns are both fine as is
    # (have double-checked that node ids are all automatically and correctly
    # updated by the tskit.TableCollection.union() method)

    # metadata column is currently unused by both msprime and gnx


    # reformat sites table content:
    #-----------------------------

    # double-check that position and ancestral state columns are identical in
    # the recipient Species' TableCollection;
    # if so, then just duplicate that table here
    # (because the only difference between the two in the metadata column,
    # which in gnx stores whether a site started as a neutral ('n')
    # or non-neutral ('t', for trait-influencing) site
    # NOTE: not running these checks for msprime-derived tables because
    #       they can lack some sites, depending on recomb rate, Ne, etc
    if source_software == 'gnx':
        for col in ['position', 'ancestral_state']:
            assert np.all(getattr(tc.sites, col) ==
                          getattr(recip_spp._tc.sites, col)), (
                                f"The {col} column in the sites table "
                                f"generated by msprime disagrees with the {col} "
                                 "column in the recipient Species' sites table.")
    tc.sites.replace_with(recip_spp._tc.sites)


    # reformat mutations table content:
    #---------------------------------

    # site, node, and parent columns are the core pieces needed to recover the
    # genotypes produced by the mutations, and we want whatever genotypes
    # msprime has simulated, so these columns should all remain untouched
    # site vals are all in interval [0, L-1]
    assert (np.all(0 <= tc.mutations.site) and
            np.all(tc.mutations.site <= (recip_spp.gen_arch.L-1))), ("The "
                        "site column in the mutations table contains values "
                        "not in the interval [0, recip_spp.gen_arch.L-1].")

    # Times for mutations created at the start of a gnx simulation, to
    # match the starting allele frequency array ('p') provided to a gnx
    # model, are artificial as they are all forcibly mutated at at once.
    # Thus, those times are NaN. Unfortunately, tskit does not all
    # mutation times to be mixed known and unknwon, so that means we must
    # set all mutation times to NaN
    time = [tskit.UNKNOWN_TIME] * len(tc.mutations.time)
    tc.mutations.time = time

    # double-check that the derived_state_offset column is just serial
    # integers from 0 to the table length, inclusive
    assert np.all(tc.mutations.derived_state_offset ==
                  np.linspace(0,
                              tc.mutations.num_rows,
                              tc.mutations.num_rows + 1)), ("The "
                "derived_state_offset column of the mutations table should "
                "contain serial integers from 0 to the table length, "
                "inclusive, but does not.")

    # metadata, metadata_offset, and metadata_schema are all currently unused
    # in both msprime and gnx


    # reformat migrations table content:
    #-----------------------------------

    # this table is not used in gnx, but if it was used in a complex msprime
    # Demography that simulated individuals for introduction into gnx then it
    # will not be empty in this TableCollection; in that case, its contents
    # will just be unioned into the otherwise empty migrations table of the 
    # recipient Species and then left as is


    # reformat populations table content:
    #-----------------------------------

    # This table is unused by gnx, but the msprime table has
    # basic population names (e.g., 'pop_0') and descriptions ('' by default,
    # but a user could specify).
    # We won't touch or do anything with this data, and the
    # TableCollection.union() method will pull it along without editing,
    # but the rows in this table will be appended to existing rows in the gnx
    # model and the nodes table will be constructed so as to match those indices

    # reformat provenances table content:
    #-----------------------------------

    # no need to do anything here; all provenances info will be retained and
    # combined when the union() method is called on the recipient Species'
    # TableCollection

    # return the gnx idxs that were set in the individuals table, for
    # downstream use
    return gnx_individ_idxs


def _get_lineage_dicts(spp, nodes, loci, t_curr, drop_before_sim=True,
                       time_before_present=True, max_time_ago=None,
                       min_time_ago=None):
    '''
    For a sample of nodes and a sample of loci, get the
    loci-sequentially ordered dict of lineage dicts
    (containing each node in a child node's lineage as the
    key and its birth-time and birth-location as a length-2 list of values

    Thus, the returned value is of the structure:
        {locus_id_1: {node_id_1: {curr_node_id: (birth_time,
                                             array([birth_loc_par_1,
                                                    birth_loc_par_2])),
                              prev_node_id: (birth_time,
                                             array([birth_loc_par_1,
                                                    birth_loc_par_2])),
                              .
                              .
                              .
                              oldest_node_id: (birth_time,
                                             array([birth_loc_par_1,
                                                    birth_loc_par_2]))},
                    .
                    .
                    .
                    node_id_last: {curr_node_id: ...
                                    ...
                                    }},
        .
        .
        .
        locus_id_last: {node_id_1: {curr_node_id: ...
                                    ...
                                    }}
        }
    '''

    # make sure the TableCollection is sorted, then get the TreeSequence
    try:
        ts = spp._tc.tree_sequence()
    except Exception:
        raise Exception(("The species' TableCollection must be sorted and "
                         "simplified before this method is called."))
    # get the tree numbers corresponding to each locus
    treenums = _get_treenums(ts, loci)
    # get the trees
    trees = ts.aslist()
    # create the output data structure
    lineage_dicts = {}
    # get the lineage_dict for each locus
    for loc, treenum in zip(loci, treenums):
        # store the lineage dicts for the current locus
        lineage_dicts_curr_loc = _get_lineage_dicts_one_tree(spp._tc,
                                    trees[treenum], nodes, t_curr,
                                    drop_before_sim=drop_before_sim,
                                    time_before_present=time_before_present,
                                    max_time_ago=max_time_ago,
                                    min_time_ago=min_time_ago)
        lineage_dicts[loc] = lineage_dicts_curr_loc
    return lineage_dicts


def _get_lineage_dicts_one_tree(tc, tree, nodes, t_curr, drop_before_sim=True,
                                time_before_present=True, max_time_ago=None,
                                min_time_ago=None):
    '''
    For a sample of nodes, and for a given tree,
    get dicts of each node's lineage nodes with their birth locs and birth
    times
    '''

    #create an output master dict
    lineage_dicts = {}
    # get each sample node's lineage dict
    for node in nodes:
        lineage = _get_lineage(tree, [node])
        lineage_dict =dict(zip(lineage, _get_lineage_times_and_locs(
                                    tc, lineage, t_curr,
                                    drop_before_sim=drop_before_sim,
                                    time_before_present=time_before_present)))

         # filter for time, if requested
        if max_time_ago is not None or min_time_ago is not None:
            # if only one time limit was set, set the other correctly, based on
            # whether or not time is expressed as time before present
            if max_time_ago is None:
                max_time_ago = np.inf
            if min_time_ago is None:
                min_time_ago = -np.inf
            # filter by time
            lineage_dict = {node: val for node, val in lineage_dict.items() if (
                            min_time_ago <= val[0] <= max_time_ago)}

        lineage_dicts[node] = lineage_dict

    return lineage_dicts


# get birth-locations and times of nodes along a pedigree
def _get_lineage_times_and_locs(tc, lineage, t_curr, drop_before_sim=True,
                                time_before_present=True):
    locs = []
    times = []
    for node in lineage:
        time = tc.nodes[node].time
        # NOTE: optionally drop all nodes from before the simulation
        # (i.e. nodes which were added to TableCollection by
        # backward-time simulation with ms); THIS IS THE DEFAULT
        # (or include all, if drop_before_sim is False)
        if (drop_before_sim and time < 0) or not drop_before_sim:
            times.append(time)
            individ = tc.nodes[node].individual
            loc = tc.individuals[individ].location
            locs.append(loc)
    # express times as time before present, if requested
    if time_before_present:
        # NOTE: need to make t_curr negative, to reflect the fact that time is
        # stored as positive integeres in Geonomics data structures but as
        # negative integers in tskit TableCollection (because they require time
        # expressed as time before present, so using positive integers would
        # not allow models to be run indefinitely, for exploratory/fun purposes
        times = [t - -t_curr for t in times]
    return(zip(times, locs))


def _get_lineage(tree, nodes):
    '''
    Recursive algorithm for getting all parents of a node
    '''
    parent = tree.parent(nodes[-1])
    if parent == tskit.NULL:
        return(nodes)
    else:
        nodes.append(parent)
        return(_get_lineage(tree, nodes))


# get a dict of the interval (i.e. tree) number
# keyed to each locus in a list of loci
def _get_treenums(ts, loci, as_dict=False):
    interval_right_ends = [t.get_interval()[1] for t in ts.trees()]
    treenums = [bisect.bisect_left(interval_right_ends, loc) for loc in loci]
    if as_dict:
        treenums = {loc: treenum for loc, treenum in zip(loci, treenums)}
    return(treenums)


# calculate a stat for the lineage corresponding to the provided lineage_dict
def _calc_lineage_stat(lin_dict, stat):
    assert stat in ['dir', 'dist', 'time', 'speed'], ("The only valid "
                                                       "statistics are: "
                                                       "direction ('dir'), "
                                                       "distance ('dist'), "
                                                       "time ('time'), "
                                                       "and speed ('speed').")
    fn_dict = {'dir': _calc_lineage_direction,
               'dist': _calc_lineage_distance,
               'time': _calc_lineage_time,
               'speed': _calc_lineage_speed
              }
    val = fn_dict[stat](lin_dict)
    return val


# calculate the angular direction of gene flow along a locus' lineage
def _calc_lineage_direction(lin_dict):
    """
    returns 'compass' directions of gene flow, with directions expressed
    as angles (in degrees) clockwise from compass north
    """
    # if the lineage has 1 or less node, don't calculate anything
    if len(lin_dict.values()) < 2:
        return
    else:
        # get x and y distances between beginning and ending points
        # (begins at the bottom of the lineage dict, furthest back in time)
        beg_loc = [*lin_dict.values()][-1][1]
        # (ends at the top of the lineage dict, in the current time step)
        end_loc = [*lin_dict.values()][0][1]
        # get the difference in the x and y axes between the lineage's
        # beginning location and its ending (i.e. current) location
        # NOTE: using range(2) will grab only the first location values,
        # which is perfect, because models with traits also save phenotypic
        # and fitness values in the individuals table's 'location' column
        x_diff, y_diff = [end_loc[i] - beg_loc[i] for i in range(2)]
        # get the counterclockwise angle, expressed in degrees
        # from the vector (X,Y) = (1,0),
        # with 0 to 180 in quadrants 1 & 2, 0 to -180 in quadrants 3 & 4
        ang = np.rad2deg(np.arctan2(y_diff, x_diff))
        # convert to all positive values, 0 - 360
        if ang < 0:
            ang += 360
        #convert to all positive clockwise angles,
        # expressed as 0 at compass north
        # (i.e. angles clockwise from vector (X,Y) = (0,1)
        ang = (-ang + 90) % 360
        return ang


# calculate the geographic distance (in cell widths)
# of gene flow along a locus' lineage
def _calc_lineage_distance(lin_dict):
    if len(lin_dict.values()) < 2:
        return
    else:
        # get x and y distances between beginning and ending points
        beg_loc = [*lin_dict.values()][0][1]
        end_loc = [*lin_dict.values()][-1][1]
        x_diff, y_diff = [end_loc[i] - beg_loc[i] for i in range(2)]
        #Pythagoras
        dist = np.sqrt(x_diff**2 + y_diff**2)
        return dist


# calculate the total time to the simulation's MRCA of a locus' lineage
def _calc_lineage_time(lin_dict):
    if len(lin_dict.values()) < 2:
        return
    else:
        beg_time = [*lin_dict.values()][0][0]
        end_time = [*lin_dict.values()][-1][0]
        time = end_time - beg_time
        return time


# calculate the speed of gene flow in a lineage (in cell widths/time steps)
def _calc_lineage_speed(lin_dict):
    if len(lin_dict.values()) < 2:
        return
    else:
        dist = _calc_lineage_distance(lin_dict)
        time = _calc_lineage_time(lin_dict)
        speed = dist/time
        return speed




#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################


import numpy as np
import geonomics as gnx
import msprime
import tskit
import warnings
from copy import deepcopy
from collections import Counter as C

