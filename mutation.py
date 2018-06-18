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


def estimate_total_num_mutations(params, pop):
    #NOTE: this would actually be a pretty poor estimate, because mutations will occur in new individuals, not some static population

    mean_births = np.mean(pop.n_births[-params['model']['burn']['burn_T_min']:])
    est = mean_births * pop.genomic_arch.L * params['model']['main']['T'] * pop.genomic_arch.mu

    #give a decent overestimate

    est = int(2.5 * est)

    return(est)


def mutate(pop, log = False):

    newborns = {i:v.age for i,v in pop.individs.items() if v.age == 0}

    mutation = r.binomial(1, pop.genomic_arch.mu*pop.genomic_arch.L*len(newborns))


    if mutation == 1:
        #randomly choose an individual
        ind = r.choice(list(newborns.keys()))
        #randomly choose a locus from among the mutables
        r.shuffle(pop.genomic_arch.mutables)
        loc = pop.genomic_arch.mutables.pop()
        #randomly choose a chromosome
        chrom = r.binomial(1,0.5)

        #then mutate this individual at this locus
        pop.individs[ind].genome[loc,chrom] = 1




        #NOTE: Change this to something more generalizable in the main script
        if log == True:
            with open('./mutation_log.txt', 'a') as f:
                f.write('locus: %i,t: %i' % (locus, t))
        
        print('#\n#\n#\n#\n#\n#\n#\n')
        print('MUTATION:\n\t INDIVIDUAL %i,  LOCUS %i\n\t timestep %i\n\n' % (ind, loc, pop.t))
        print('#\n#\n#\n#\n#\n#\n#\n')

