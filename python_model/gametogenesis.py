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



def segregate():
        segregation = r.binomial(1, 0.5) #determine which homologous chromosome will segregate to the gamete (i.e. whether to start the recombination process on the left or right col of the genome array)

        return segregation


#function to generate a gamete from given individual.Individual and genome.Genomic_Architecture instances
def gametogenerate(individ, genomic_arch):
    gamete = {}
    for i, c in individ.genome.genome.items():
        recombinant = recombination.recombine(chromosome_num = i, genomic_arch = genomic_arch)
        segregation = segregate()
        chromatid_path = np.array([int(segregation!=locus) for locus in recombinant]) #set path to match whichever of the two homologous chromosomes should be taken, per the variable 'segregation' above
        chromatid = np.array([c[i_j] for i_j in zip(range(len(c)),chromatid_path)])  #produce the gamete resulting from the chromatid 'path'
        gamete[i] = chromatid

    return gamete


