#!/usr/bin/python
#mutation.py

'''
##########################################

Module name:                  mutation

Module contains:
                              - function for simulating mutation across the genome, according to input parameters
                              - associated functions


Author:                       Drew Ellison Hart
Email:                        drew.hart@berkeley.edu
Github:                       URL
Start date:                   12-28-15
Documentation:                URL


##########################################
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

def calc_estimated_total_mutations(params, pop):
    #NOTE: this would actually be a pretty poor estimate, because mutations will occur in new individuals, not some static population
    mean_births = np.mean(pop.n_births[-params.model.its.burn.T_min:])
    est = mean_births * pop.gen_arch.L * params.model.its.main.T * pop.gen_arch.mu_tot
    #give a decent overestimate
    est = int(2.5 * est)
    return est 


def do_neutral_mutation(pop, offspring, locus=None, individ=None):
    #randomly choose an individual to mutate, if not provided
    if individ is None:
        individ = r.choice(offspring)
    else:
        assert individ in pop.keys(), 'ERROR: The individual provided is not in pop.keys().'
    #randomly choose a locus from among the mutables, if not provided
    if locus is None:
        locus = r.choice([*pop.gen_arch.mutable_loci])
    else:
        assert locus in pop.gen_arch.mutable_loci, 'ERROR: The locus provided is not in the pop.gen_arch.mutable_loci set.'
    #remove the locus from the mutable_loci set
    pop.gen_arch.mutable_loci.remove(locus)
    #randomly choose a chromosome on which to place the mutation in the individual
    chrom = r.binomial(1,0.5)
    #then mutate this individual at this locus
    pop[individ].genome[locus,chrom] = 1
    return(individ, locus)


def do_nonneutral_mutation(pop, offspring, locus=None, individ=None):
    #choose a new locus, if not provided
    if locus is None:
        locus = pop.gen_arch.draw_mut_loci()
        assert locus in pop.gen_arch.mutable_loci, 'ERROR: The locus provided is not in pop.gen_arch.mutable_loci.'
    #remove the locus from the mutable_loci and neut_loci sets
    pop.gen_arch.mutable_loci.remove(locus)
    pop.gen_arch.neut_loci.remove(locus)
    #add the locus to the nonneut_loci set
    pop.gen_arch.nonneut_loci.update({locus})
    #draw an individual to mutate from the list of offspring provided, unless individ is provided
    if individ is None:
        individ = r.choice(offspring)
    #TODO: Need to check that the individual provided is among the list of offspring? Or make offspring an optional arg too?
    #create the mutation
    pop[individ].genome[locus,r.binomial(1, 0.5)] = 1 
    #return the locus and individual
    return(locus, individ)


def do_trait_mutation(pop, offspring, trait_num, alpha=None, locus=None, individ=None):
    #run the do_nonneutral_mutation function, to select locus and individ (unless already provided) and update
    #the mutable_loci, neut_loci, and nonneut_loci sets, and change one of this locus' alleles to 1 in the mutated individual
    locus, individ = do_nonneutral_mutation(pop = pop, offspring = offspring, locus = locus, individ = individ)
    #choose an effect size for this locus, if not provided
    if alpha is None:
        alpha = pop.gen_arch.draw_trait_alpha(trait_num)
    #update the trait's loci and alpha attributes accordingly
    pop.gen_arch.set_trait_loci(trait_num, mutational = True, loci = locus, alpha = alpha)
    #return the locus and individual
    return(locus, individ)


def do_deleterious_mutation(pop, offspring, locus=None, s=None, individ=None):
    #run the do_nonneutral_mutation function, to select locus and individ (unless already provided) and update
    #the mutable_loci, neut_loci, and nonneut_loci sets, and change one of this locus' alleles to 1 in the mutated individual
    locus, individ = do_nonneutral_mutation(pop = pop, offspring = offspring, locus = locus, individ = individ)
    #choose a selection coefficient for this locus, if not provided
    if s is None:
        s = pop.gen_arch.draw_delet_s()
    #update the pop.gen_arch.delet_loci OrderedDict
    pop.gen_arch.delet_loci.update({locus: s})
    #return the locus and individual
    return(locus, individ)


def do_planned_mutation(planned_mut_params):
    pass


def do_mutation(offspring, pop, log=None):
    #draw number of mutations from a binomial trial with number of trials equal to
    #len(offspring)*pop.gen_arch.L and prob = sum(all mutation rates)
    n_muts = r.binomial(n = len(offspring) * pop.gen_arch.L, p = pop.gen_arch.mu_tot)
    print('# MUTS : %i' % n_muts)

    #if n_muts >0, draw which type of mutation each resulting mutation should be 
    #(weighted by relative mutation rates)
    if n_muts > 0:
        muts = pop.gen_arch.draw_mut_types(num = n_muts)

        #for each mutation, add the corresponding mutation function to a queue
        mut_queue = []
        for mut in muts:
            #add the appropriate mut_fn to the queue
            mut_fn = pop.gen_arch.mut_fns[mut]
            mut_queue.append(mut_fn)

        #add to the queue any planned mutations, if necessary
        if pop.gen_arch.planned_muts:
            if pop.gen_arch.planned_muts.next[0] <= pop.t:
                mut_queue.append(planned_muts.next[1])
                pop.gen_arch.planned_muts.set_next()

        #execute all functions in the queue
        for fn in mut_queue:
            individ, locus = fn(pop, offspring)
            log_msg = 'MUTATION:\n\t INDIVIDUAL %i,  LOCUS %i\n\t timestep %i\n\n' % (individ, locus, pop.t)
            if log:
                with open(log, 'a') as f:
                    f.write(log_msg)
            print(log_msg)

