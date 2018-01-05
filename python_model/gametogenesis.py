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

#function to give a recombination path (i.e. list of homologous chromosome IDs, as 0s and 1s, for all loci, from which 
#a new gamete's genotype is to be drawn) for a given chromosome number and genomic architecture's genetic map for that chromosome
def recombine(genomic_arch):
    recombination = r.binomial(1, genomic_arch.r) #determine all recombination events (i.e. 1's)
    recombination = np.cumsum(recombination)%2 #determine a recombination 'path' down the genome array (i.e. which col to pull from each row)

    return recombination



#function to generate a gamete from given individual.Individual and genome.Genomic_Architecture instances
def gametogenerate(individ, genomic_arch):
    recombination_path = recombine(genomic_arch = genomic_arch)
    #produce the gamete resulting from the recombination path
    gamete = individ.genome[range(genomic_arch.L),recombination_path]
    #gamete = np.array([individ.genome[n][recombination_path[n]] for n in range(len(recombination_path))])  
    return gamete


