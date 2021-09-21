#!/usr/bin/python
# mating.py

# flake8: noqa

'''
Functions to implement mating operations.
'''

#other imports
from scipy.spatial import cKDTree
import numpy as np
import numpy.random as r
from operator import itemgetter as ig
from itertools import repeat, starmap


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

def _find_mates(spp, sex=False, repro_age=None, choose_nearest=False,
                inverse_dist_mating=False):
    # NOTE: In an IDEAL world, this would work as follows:
    # for all individuals, find set of nearest neighbors (excluding selves)
    # for females:
        # choose male from among those nn's, based on:
            # probabilistic function of distance
        # if sexual_selection:
            # (product of?) selection coefficients
    ######################################################
    # First, query the scipy.spatial.cKDTree (i.e. spp._kd_tree) for
    #nearest-neigh pairs
    pairs = spp._get_mating_pairs(choose_nearest=choose_nearest,
                                  inverse_dist_mating=inverse_dist_mating)
    ####################################################
    # Then, operationalize sexes, if being used, and find all 
    #available pairs within max distance
    if sex:
        # np.array of the sexes of all individuals
        sexes = np.array([ind.sex for ind in spp.values()])
        # array of couplings for all females with 
        #nearest individual < mating_radius
        # i.e.  [AT LEAST 1 INDIVID < mating_radius (OTHERWISE 
        #scipy.spatial.cKDTree associates them with a nonexistent 
        #neighbour and infinite distance if choose_nearest is True,
        #otherwise they're associated with themselves at distance==0)]
        #& [female]
        available_females =np.array(sexes[pairs[:,0]] == 0)
        # check that nearest individuals to paired females are males,
        #otherwise remove
        mating_pairs = pairs[available_females][sexes[pairs[:,
                                            1][available_females]] == 1]
    else:
        # delete inverse-equal pairs from list (i.e. if 3-paired-to-5 and 
        #5-paired-to-3 are both listed, drop one)
        # NOTE: only if mating_radius is not None, because if it is None then
        #       WF-style panmixia should be enforced, in which case we should
        #       allow for draws with replacement to generate duplicate pairs
        if spp.mating_radius is not None:
            mating_pairs = np.array(list(map(tuple, set(map(frozenset, pairs)))))
        else:
            mating_pairs = pairs
    ##############################################
    # Check if any available pairs have been found thus far, proceed if so,
    #otherwise return empty array of pairs
    if len(mating_pairs) == 0:
        mates = np.array([])
        return mates
    else:
        pass
    ###############################################
    # Then, operationalize age-structuing, if being used
    # NOTE: I have NOT restricted this so that each individual can only 
    #mate once per generation, but I should probably figure that out and 
    #add it here!
    if (repro_age is not None
        and np.any(np.atleast_1d(repro_age) > 0)):
        # np.array of the ages of all individuals
        ages = np.array([ind.age for ind in spp.values()])
        # if sexual species, repro_age expected to be a tuple or list of 
        #numerics of length 2
        if sex:
            assert repro_age.__class__.__name__ in ('ndarray', 'list',
                'tuple'), ("For a sexual and age-structured species, "
                "age at sexual maturity, 'repro_age', must be expressed as an "
                "iterable of numerics (i.e. a list, tuple, or numpy.ndarray "
                "of floats or integers); you have provided a:\n\t"
                "%s") % repro_age.__class__.__name__
            yes_f = np.array(ages[mating_pairs[:, 0]] >= repro_age[0])
            yes_m = np.array(ages[mating_pairs[:, 1]] >= repro_age[1])
            mating_pairs = mating_pairs[yes_f & yes_m]
        else:  # if not sexual species, repro_age expected to be a numeric
            assert repro_age.__class__.__name__ in ('int',
                'float'), ("For a non-sexual and age-structured species, "
                "the age at sexual maurity, 'repro_age', must be expressed as "
                "a single, non-iterable numeric (i.e.  float or integer); you "
                "have provided a:\n\t%s") % repro_age.__class__.__name__
            yes_0 = np.array(ages[mating_pairs[:, 0]] >= repro_age)
            yes_1 = np.array(ages[mating_pairs[:, 1]] >= repro_age)
            yes = np.sum(ages[mating_pairs] >= repro_age, axis = 1) == 2
            mating_pairs = mating_pairs[yes]
    #########################################
    # Then, if at least one valid pair was found, match pairs' indices
    # (in the data structure returned from get_mating_pairs)
    # back to their individual ids
    if len(mating_pairs) > 0:
        # finally, link individuals' ordinal indices back to the initially
        #created structure, to get individuals' proper keys
        f = ig(*mating_pairs.flatten())
        mates = np.array(f([*spp])).reshape(mating_pairs.shape)
    else:
        mates = np.array([])
    # Return an array or arrays, each inner array containing a mating-pair
    return mates


def _draw_n_births(num_pairs, n_births_distr_lambda):
    #NOTE: subtracting nearly 1 from the lambda, then adding 1 to the drawn
    #value, guarantees that at least 1 offspring will be born for each pair
    #chosen to mate
    num_births = r.poisson((n_births_distr_lambda), num_pairs)
    num_births = np.clip(num_births, a_min = 1, a_max = None)
    return num_births


# function for mating for a single offspring for a certain mating pair
def _do_mating_sngl_offspr(spp, pair, recomb_keys):
    # choose a start homologue
    # NOTE: for now, this will only work for diploidy
    start_homologues = r.binomial(1, 0.5, 2)

    # if using tskit, generate the segment info for each parent's new,
    # recombined gamete
    # NOTE: reverse the order of the recomb_paths object so that the resulting
    # segment info objects' order matches the order in which the recombination
    # paths are popped in the following line of code (because they will pop
    # off the right end of their list)
    if spp.gen_arch.use_tskit:
        seg_info = [spp.gen_arch.recombinations._get_seg_info(
                        start_homologue=start_homologues[i],
                        event_key=k,
                        node_ids=np.array([spp[pair[i]]._nodes_tab_ids[
                             homol] for homol in range(
                             spp.gen_arch.x)])
                        ) for i, k in enumerate(recomb_keys)]
    else:
        seg_info = None

    # generate each parent's new, recombined gamete (and make into a new genome
    # by stacking and transposing)
    # NOTE: the first path winds up the left side of the new genome, thus
    # homologue 0; this ensures that the first segment-info object in the
    # `seg_info` list and the first (i.e. left) half of the new genome
    # both correspond to the genomic material inherited from the first of the
    # two parents whose ids are listed in `pair`
    if ((spp.gen_arch.use_tskit and len(spp.gen_arch.nonneut_loci) > 0)
        or not spp.gen_arch.use_tskit):
        subsetters = [spp.gen_arch.recombinations._get_subsetter(
                                        event_key=k) for k in recomb_keys]
        #NOTE: flip the genome L-R before subsetting, if the start homologue is
        # 1, then flatten genome and subset
        new_genome = [np.fliplr(spp[ind].g).flatten(
            )[np.bool8([*sub])] if hom else spp[ind].g.flatten(
            )[np.bool8([*sub])] for ind, hom, sub in zip(pair,
                                               start_homologues, subsetters)]
        new_genome = np.vstack(new_genome).T
    else:
        new_genome = None
    return new_genome, seg_info


# function for mating a chosen mating-pair
def _do_mating_sngl_pair(spp, pair, n_births, pairs_recomb_keys):
    # generate a list of n offspring for the given pair
    offspring = [_do_mating_sngl_offspr(spp, pair,
                                        [pairs_recomb_keys.pop(
                                        ) for _ in range(2)]) for off in range(
                                                                  n_births)]
    return offspring


# function for mating a whole species
def _do_mating(spp, mating_pairs, n_births, recomb_keys):
    """
    NOT INTENDED FOR PUBLIC USE

    spp: a gnx Species object

    mating_pairs: an nx2 2d np.array of ints, where n = number of mating pairs,
                  2 columns correspond to each of the two mates in each pair,
                  and each int is the index number of an Individual in a pair

    n_births: a length-n 1d np.array of ints, where n = number of mating pairs
              and each int is the number of offspring the nth pair will produce

    recomb_keys: an m-length list of ints, where m = twice the
                 total number of births and each int is a randomly drawn
                 recombination event's key
    """
    # use list of number of offspring per mating pair to create
    # list of start and stop indices for subsetting the recomb paths into
    # groups of paths to be used for each gamete-production event for each pair
    start_stop_idxs = np.hstack((0, np.cumsum([2*n for n in n_births])))
    # get the recombination paths to be used by each pair
    pairs_recomb_keys = [recomb_keys[start_stop_idxs[i]: start_stop_idxs[
                                    i + 1]] for i in range(len(n_births))]
    # use the pairs' paths to produce recombinant genomes for each pair (where
    # number of genomes equals that pair's number of offspring)
    new_genomes = list(starmap(_do_mating_sngl_pair, zip(repeat(spp),
                               mating_pairs, n_births, pairs_recomb_keys)))
    return new_genomes
