#!/usr/bin/python
#gametogenesis.py

'''
##########################################

Module name:              gametogenerate

Module contains:
                          - function for generating a gamete
                          - associated functions


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               12-28-15
Documentation:            URL


##########################################
'''

import numpy as np
import numpy.random as r





#------------------------------------
# CLASSES ---------------------------
#------------------------------------





#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------


#function to generate a gamete from given individual.Individual and genome.Genomic_Architecture instances
def gametogenerate(individ, recomb_path):
    #produce the gamete resulting from the recombination path
    gamete = individ.genome.flatten()[recomb_path]
    return(gamete)



#carry out recombination, using the lookup array in a Genomic_Architecture object
def recombine(r_lookup, n_recombinants):
    recombinants = np.array([r.choice(r_lookup[i,], size = n_recombinants, replace = True) for i in range(len(r_lookup))])
    recombinants = np.cumsum(recombinants, axis = 0)%2
    return(recombinants)


