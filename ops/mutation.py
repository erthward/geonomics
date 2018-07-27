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


def calc_estimated_total_mutations(params, pop):
    #NOTE: this would actually be a pretty poor estimate, because mutations will occur in new individuals, not some static population

    mean_births = np.mean(pop.n_births[-params.model.its.burn.T_min:])
    est = mean_births * pop.gen_arch.L * params.model.its.main.T * pop.gen_arch.mu

    #give a decent overestimate

    est = int(2.5 * est)

    return(est)


def do_mutation(pop, log = False):

    newborns = {i:v.age for i,v in pop.items() if v.age == 0}

    mutation = r.binomial(1, pop.gen_arch.mu*pop.gen_arch.L*len(newborns))


    if mutation == 1:
        #randomly choose an individual
        ind = r.choice(list(newborns.keys()))
        #randomly choose a locus from among the mutables
        r.shuffle(pop.gen_arch.mutables)
        loc = pop.gen_arch.mutables.pop()
        #randomly choose a chromosome
        chrom = r.binomial(1,0.5)

        #then mutate this individual at this locus
        pop[ind].genome[loc,chrom] = 1




        #NOTE: Change this to something more generalizable in the main script
        if log == True:
            with open('./mutation_log.txt', 'a') as f:
                f.write('locus: %i,t: %i' % (locus, t))
        
        print('#\n#\n#\n#\n#\n#\n#\n')
        print('MUTATION:\n\t INDIVIDUAL %i,  LOCUS %i\n\t timestep %i\n\n' % (ind, loc, pop.t))
        print('#\n#\n#\n#\n#\n#\n#\n')

