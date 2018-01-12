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


#####function to give a recombination path (i.e. list of homologous chromosome IDs, as 0s and 1s, for all loci, from which 
####a new gamete's genotype is to be drawn) for a given chromosome number and genomic architecture's genetic map for that chromosome
###def recombine(genomic_arch):
###    recombination = r.binomial(1, genomic_arch.r) #determine all recombination events (i.e. 1's)
###    recombination = np.cumsum(recombination)%2 #determine a recombination 'path' down the genome array (i.e. which col to pull from each row)
###
###    return recombination



#function to generate a gamete from given individual.Individual and genome.Genomic_Architecture instances
def gametogenerate(individ, recomb_path):
    #produce the gamete resulting from the recombination path
    genome = individ.genome
    gamete = genome[range(np.shape(genome)[0]),recomb_path]
    #gamete = np.array([individ.genome[n][recombination_path[n]] for n in range(len(recombination_path))])  
    return gamete



#carry out recombination, using the lookup array in a Genomic_Architecture object
def recombine(la, n_recombinants):
    recombinants = np.array([r.choice(la[i,], size = n_recombinants, replace = True) for i in range(len(la))])
    recombinants = np.cumsum(recombinants, axis = 0)%2
    return(recombinants)


