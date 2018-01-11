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


import numpy as np
from numpy import random as r
import random
  

#------------------------------------
# CLASSES ---------------------------
#------------------------------------





#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------


def est_total_num_mutations(params, pop):
    est_total_num_muts = int(params['mu'] * params['L'] * sum(pop.K.raster) * params['T'])
    return(est_total_num_muts)


def mutate(pop, individ, genomic_arch, t, alpha_mut_s = 25, beta_mut_s = 0.5):

    #NOTE: should I incorporate the number of offspring produced in the following calculation? or should I make this a method of individuals, and do it once for every offspring?

    mutation = r.binomial(1, genomic_arch.mu*genomic_arch.L)


    if mutation:
        #randomly choose a locus from among the neutral loci
        locus = random.choice(np.array(range(genomic_arch.L))[genomic_arch.non_neutral])



        #NOTE: Change this to something more generalizable in the main script
        with open('./mutation_log.txt', 'a') as f:
            f.write('locus: %i,t: %i' % (locus, t))
            
            
        
        print 'MUTATION AT: locus %i' % (locus)
        raw_input()



        #set all of the current population's alleles at this locus to 0
        #NOTE: This is a huge cheat! But for now seems to most obvious way to implement this
            #DEH 12-04-18: and shouldn't matter much because I am going to overhaul the whole mutation approach and module
        for ind in pop.individs.values():
            ind.genome[locus,:] = [0,0]

        #insert mutant in offspring's genome
        mutant = [0,1]
        r.shuffle(mutant)
        individ.genome[locus,:] = mutant

        ###NOTE DEH 12-04-18: This is now completely obsolete
        ###set selection coefficient to (probably, b/c drawn from right-heavy beta distribution) something highly advantageous
        ###genomic_arch.s[chrom][locus] = r.beta(alpha_mut_s, beta_mut_s)

        ###set locus to be non-neutral
        ###genomic_arch.non_neutral[chrom][locus] = True
        
        





#NOTE: use the above to write a function that INTENTIONALLY introduces a mutation
   
def introduce_mutation(pop):
    pass



