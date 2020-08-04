#!/usr/bin/python
#mutation.py

'''
Functions to implement mutation operations
'''

#other imports
import numpy as np
from numpy import random as r
import random
import re


#------------------------------------
# CLASSES ---------------------------
#------------------------------------


#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------

def _calc_estimated_total_mutations(spp, burn_T, T):
    #NOTE: this would actually be a pretty poor estimate, because mutations
    #will occur in new individuals, not some static population
    mean_births = np.mean(spp.n_births[-burn_T:])
    est = mean_births * spp.gen_arch.L * T * spp.gen_arch._mu_tot
    #give a decent overestimate
    est = int(2.5 * est)
    return est

# add a row to the tskit.TableCollection.mutations table for a new mutation
def _do_add_row_muts_table(spp, individ, homol, locus):
    # update the tskit.TableCollection.mutations table
    node_id = spp[individ]._nodes_tab_ids[homol]
    mut_id = spp._tc.mutations.add_row(site=locus, node=node_id,
                                       derived_state='1')
    return mut_id


# do a neutral mutation for a single Individual, chosen from the offspring
def _do_neutral_mutation(spp, offspring, locus=None, individ=None):
    #randomly choose an individual to mutate, if not provided
    if individ is None:
        individ = r.choice(offspring)
    else:
        assert individ in spp.keys(), ('The individual provided '
                                       'is not in spp.keys().')
    #randomly choose a locus from among the mutables, if not provided
    if locus is None:
        locus = spp.gen_arch._mutables.pop()
    else:
        assert locus in spp.gen_arch._mutables, ('The locus provided '
                                                 'is not in the list of '
                                                 'valid mutable loci '
                                                 '(spp.gen_arch._mutables).')
    #randomly choose a homologue on which to place the mutation
    #in the individual
    homol = r.binomial(1,0.5)
    # NOTE: not mutating the individual's genome because we're only tracking
    # non-neutral mutations now

    # add a row to the tskit.TableCollection.mutations table
    _do_add_row_muts_table(spp, individ, homol, locus)
    return(individ, locus)


# do a non-neutral mutation for a single Individual, chosen from the offspring
def _do_nonneutral_mutation(spp, offspring, locus=None, individ=None,
                           trait_nums=None, delet_s=None):
    # choose a new locus, if not provided
    if locus is None:
        locus = spp.gen_arch._mutables.pop()
    # otherwise, remove from the mutables the locus that was provided
    else:
        spp.gen_arch._mutables.remove(locus)
    # draw an individual to mutate from the list of offspring provided,
    # unless individ is provided
    if individ is None:
        individ = r.choice(offspring)
    else:
        assert individ in offspring, ('Individual provided for mutation (%i) '
                                      'does not appear to be among the '
                                      'list of current-timestep '
                                      'offspring (%s)') % (individ,
                                                           str(offspring))
    # add the locus to the gen_arch.nonneut_loci array
    idx = spp.gen_arch._add_nonneut_locus(locus, trait_nums, delet_s)
    #create the mutation
    homol = r.binomial(1, 0.5)
    # add a new row for this locus, with '0' genotypes,
    # to the individuals' genotype arrays
    spp._add_new_locus(idx, locus)
    # mutate the chosen individual's chosen homologue's genotype to 1
    spp[individ].g[idx, homol] = 1
    # update the individual's phenotype
    spp._set_z_individ(individ)
    # add a row to the tskit.TableCollection.mutations table
    _do_add_row_muts_table(spp, individ, homol, locus)
    # update the recombination subsetters
    spp.gen_arch.recombinations._update_subsetters(locus, idx)
    return(individ, locus)


# do a trait mutation for a single Individual, chosen from the offspring
def _do_trait_mutation(spp, offspring, trait_nums, alpha=None,
                                        locus=None, individ=None):
    #run the do_nonneutral_mutation function, to select locus
    #and individ (unless already provided) and update
    #the mutable_loci, neut_loci, and nonneut_loci, and 
    #change one of this locus' alleles to 1 in the mutated individual
    individ, locus = _do_nonneutral_mutation(spp=spp, offspring=offspring,
                                             locus=locus, individ=individ,
                                             trait_nums=trait_nums)
    return(individ, locus)


# do a deleterious mutation for a single Individual, chosen from the offspring
def _do_deleterious_mutation(spp, offspring, locus=None, s=None, individ=None):
    #choose a selection coefficient for this locus, if not provided
    if s is None:
        s = spp.gen_arch._draw_delet_s()
    #run the do_nonneutral_mutation function, to select locus and individ
    #(unless already provided) and update the mutable_loci, neut_loci, and 
    #nonneut_loci sets, and change one of this locus' alleles to 1 in
    #the mutated individual
    individ, locus = _do_nonneutral_mutation(spp=spp, offspring=offspring,
                                             locus=locus, individ=individ,
                                             delet_s=s)
    #return the locus and individual
    return(individ, locus)


#TODO: COMPLETE THIS?
def _do_planned_mutation(planned_mut_params):
    pass


# do mutations for a list of offspring
def _do_mutation(offspring, spp, log=None):
    #draw number of mutations from a binomial trial with number of trials equal
    #to len(offspring)*spp.gen_arch.L and prob = sum(all mutation rates)
    n_muts = r.binomial(n = len(offspring) * spp.gen_arch.L,
                                    p = spp.gen_arch._mu_tot)
    #TODO: ADD A PRINT STATEMENT FOR 'verbose' MODE INSTEAD
    #if n_muts > 0:
    #   print('\t**** %i mutations occurred\n\n' % n_muts)

    #if n_muts >0, draw which type of mutation
    #each resulting mutation should be (weighted by relative mutation rates)
    if n_muts > 0:
        muts = spp.gen_arch._draw_mut_types(num = n_muts)

        #for each mutation, add the corresponding mutation function to a queue
        mut_queue = []
        for mut in muts:
            #add the appropriate mut_fn to the queue
            mut_fn = spp.gen_arch._mut_fns[mut]
            mut_queue.append(mut_fn)

        #add to the queue any planned mutations, if necessary
        if spp.gen_arch._planned_muts:
            if spp.gen_arch._planned_muts.next[0] <= spp.t:
                mut_queue.append(spp.gen_arch._planned_muts.next[1])
                spp.gen_arch._planned_muts._set_next()

        #execute all functions in the queue
        for i, fn in enumerate(mut_queue):
            individ, locus = fn(spp, offspring)
            log_msg = ('MUTATION: %s\n\t INDIVIDUAL %i,  '
                'LOCUS %i\n\t timestep %i\n\n') % (muts[0], individ,
                                                   locus, spp.t)
            if log:
                with open(log, 'a') as f:
                    f.write(log_msg)
            print(log_msg)
            mut_loci_exhausted = False
