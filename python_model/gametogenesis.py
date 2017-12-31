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



import recombination




#------------------------------------
# CLASSES ---------------------------
#------------------------------------





#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------



#function to generate a gamete from given individual.Individual and genome.Genomic_Architecture instances
def gametogenerate(individ, genomic_arch):
    gamete = {}
    for i, c in individ.genome.genome.items():
        recombination_path = recombination.recombine(chromosome_num = i, genomic_arch = genomic_arch)
        chromatid = np.array([c[n][recombination_path[n]] for n in range(len(recombination_path))])  #produce the gamete resulting from the chromatid 'path'
        gamete[i] = chromatid
    return gamete


