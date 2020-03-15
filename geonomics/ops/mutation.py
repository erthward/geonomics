#!/usr/bin/python
#mutation.py

'''
Functions to implement mutation operations
'''

#other imports
import numpy as np
from numpy import random as r
import random


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
        locus = r.choice([*spp.gen_arch._mutable_loci])
    else:
        assert locus in spp.gen_arch._mutable_loci, ('The locus provided '
                            'is not in the spp.gen_arch._mutable_loci set.')
    #remove the locus from the mutable_loci set
    spp.gen_arch._mutable_loci.remove(locus)
    #randomly choose a homologue on which to place the mutation
    #in the individual
    homol = r.binomial(1,0.5)
    #then mutate this individual at this locus
    spp[individ].g[locus, homol] = 1
    # add a row to the tskit.TableCollection.mutations table
    _do_add_row_muts_table(spp, individ, homol, locus)
    return(individ, locus)


# do a non-neutral mutation for a single Individual, chosen from the offspring
def _do_nonneutral_mutation(spp, offspring, locus=None, individ=None):
    #choose a new locus, if not provided
    if locus is None:
        locus = r.choice([*spp.gen_arch._mutable_loci])
        assert locus in spp.gen_arch._mutable_loci, ('The locus provided '
                                    'is not in spp.gen_arch.mutable_loci.')
    #remove the locus from the mutable_loci and neut_loci sets
    spp.gen_arch._mutable_loci.remove(locus)
    spp.gen_arch.neut_loci.remove(locus)
    #add the locus to the nonneut_loci set
    spp.gen_arch.nonneut_loci.update({locus})
    #draw an individual to mutate from the list of offspring provided,
    #unless individ is provided
    if individ is None:
        individ = r.choice(offspring)
    #TODO: Need to check that the individual provided is among the list of
    #offspring? Or make offspring an optional arg too?
    #create the mutation
    homol = r.binomial(1, 0.5)
    spp[individ].g[locus, homol] = 1
    # add a row to the tskit.TableCollection.mutations table
    _do_add_row_muts_table(spp, individ, homol, locus)
    return(individ, locus)


# do a trait mutation for a single Individual, chosen from the offspring
def _do_trait_mutation(spp, offspring, trait_num, alpha=None,
                                        locus=None, individ=None):
    #run the do_nonneutral_mutation function, to select locus
    #and individ (unless already provided) and update
    #the mutable_loci, neut_loci, and nonneut_loci sets, and 
    #change one of this locus' alleles to 1 in the mutated individual
    individ, locus = _do_nonneutral_mutation(spp=spp, offspring=offspring,
                                             locus=locus, individ=individ)
    #choose an effect size for this locus, if not provided
    if alpha is None:
        alpha = spp.gen_arch._draw_trait_alpha(trait_num)
    #update the trait's loci and alpha attributes accordingly
    spp.gen_arch._set_trait_loci(trait_num, mutational=True, loci=locus,
                                 alpha=alpha)
    #return the locus and individual
    return(individ, locus)


# do a deleterious mutation for a single Individual, chosen from the offspring
def _do_deleterious_mutation(spp, offspring, locus=None, s=None, individ=None):
    #run the do_nonneutral_mutation function, to select locus and individ
    #(unless already provided) and update the mutable_loci, neut_loci, and 
    #nonneut_loci sets, and change one of this locus' alleles to 1 in
    #the mutated individual
    individ, locus = _do_nonneutral_mutation(spp=spp, offspring=offspring,
                                             locus=locus, individ=individ)
    #choose a selection coefficient for this locus, if not provided
    if s is None:
        s = spp.gen_arch._draw_delet_s()
    #update the spp.gen_arch.delet_loci OrderedDict
    spp.gen_arch.delet_loci.update({locus: s})
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
        for fn in mut_queue:
            individ, locus = fn(spp, offspring)
            log_msg = ('MUTATION:\n\t INDIVIDUAL %i,  '
                'LOCUS %i\n\t timestep %i\n\n') % (individ, locus, spp.t)
            if log:
                with open(log, 'a') as f:
                    f.write(log_msg)
            print(log_msg)
            mut_loci_exhausted = False
