#!/usr/bin/python
#recombination.py

'''
##########################################

Module name:              recombination  

Module contains:
                          - function for simulating recombination between homologous chromosomes 
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


def recombine(chromosome_num, genomic_arch):
    recombination = r.binomial(1, genomic_arch.r[chromosome_num]) #determine all recombination events (i.e. 1's)
    recombination = np.cumsum(recombination)%2 #determine a recombination 'path' down the genome array (i.e. which col to pull from each row)

    return recombination


