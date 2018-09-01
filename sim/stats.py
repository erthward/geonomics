#!/usr/bin/python
#stats.py


'''
##########################################

Module name:              stats

Module contents:          - definition of Stats class (i.e. structured container 
                            for stats calculated during model run)
                          - definition of functions for calculating various stats, 
                            at specified frequencies and with specified arguments, 
                            according to the contents of the params.model.stats section


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               01-01-18
Documentation:            URL


##########################################
'''

#geonomics imports


#other imports
import numpy as np
from scipy.stats.stats import pearsonr
from collections import Counter as C
import matplotlib.pyplot as plt
import time


#------------------------------------
# CLASSES ---------------------------
#------------------------------------

class Stats:
    def __init__(self, params):

        #create a dictionary to link the stats' names in the params dicts to the functions to be called by
        #them (defined in the FUNCTIONS section of this script)
        self.function_dict = {  'Nt': calc_Nt,
                                'ld':  calc_LD,
                                'het': calc_het,
                                'maf': calc_MAF,
                                'mean_fit': calc_mean_fit,
                              }

        #create a Stats.stats object, where all of the stats calculated will be stored
        self.stats = {}
        T = params.model.its.main.T
        stats_params = params.model.stats
        for stat, stat_params in stats_params.items():
            if stat_params.calc:
                self.stats[stat] = {'data': [np.nan]*int(np.floor(T/float(stat_params.freq))),
                                    'freq': stat_params.freq,
                                    #create tuple of other, stat-specific parameters, 
                                    #to later be unpacked as arguments to the appropriate stat function
                                    'other_params': dict([(k,v) for k,v in stat_params.items() if k not in ['calc', 'freq']])
                                    }

    #create a master method, to be called each timestep, which will make a list of all stats that need to be
    #calculated that timestep (based on the calculation-frequencies provided in the params dicts),
    #and then calls the functions to calculate them all and adds the results to self.stats
    def calc_stats(self, pop, t):

        calc_list = []

        for stat,params in self.stats.items():
            if t%params['freq'] == 0:
                calc_list.append(stat)

        for stat in calc_list:
            #if 'other_params' in self.stats[stat].keys():
            self.stats[stat]['data'][t] = self.function_dict[stat](pop, **self.stats[stat]['other_params'])
            #else:
                #self.stats[stat]['data'][t] = self.function_dict[stat](pop)


#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------

#TODO: either get rid of this, or just make it check the pop.Nt list instead
def calc_Nt(pop):
    Nt = len(pop)
    return(Nt)


def calc_ld(pop, plot = False):
    
    #TODO: I should also include (either as an alternative within this fn, or as separate fn) the option to calculate D'

    #TODO: I keep getting errors like the following, which could just be due to divison of small
        #floating-point numbers, but I should figure out exactly what's going on and be sure everything checks out:
                # stats.py:117: RuntimeWarning: invalid value encountered in double_scalars

    populome = np.array([ind.genome for ind in pop.values()])
    n = np.shape(populome)[0] #num individs
    x = np.shape(populome)[2] #ploidy
    N = n*x
    L = pop.genomic_arch.L
    assert L == np.shape(populome)[1], "The length of the 1th dimension of populome doesn't equal pop.genomic_arch.L"

    r2_mat = np.zeros([L]*2)-1 # -1 serves as a 'no data' value here

    for i in range(L):
        for j in range(i+1, L):
            f1_i = np.sum(populome[:,i,:], axis = None)/(N)  #calculates freq of allele 1 at locus i
            f1_j = np.sum(populome[:,j,:], axis = None)/(N)  #calculates freq of allele 1 at locus j
            f11_ij = float(np.sum(populome[:,[i,j],:].sum(axis = 1) ==2, axis = None))/(N) #calculates freq of chroms with 1_1 haplotype at loci i and j
            D_1_1 = f11_ij - (f1_i * f1_j)
            r2 = (D_1_1**2)/(f1_i*(1-f1_i)*f1_j*(1-f1_j))
            r2_mat[i,j] = r2

    if plot == True:
        fig = plt.figure()

        ax = fig.add_subplot(1,2,1)
        #plot of LD matrix

        plt.imshow(np.clip(r2_mat, a_min = 0, a_max = None), interpolation = 'nearest')

        ax = fig.add_subplot(1,2,2)

        #plot of mean linkage values
        r2_list = [r2_mat[0,1]]
        for i in range(1,L-1):
            r2_list.append(np.mean([r2_mat[i-1,i], r2_mat[i,i+1]]))
        r2_list.append(r2_mat[L-2,L-1])
        
        plt.scatter(range(L), r2_list, c = 'red', marker = 'o', s=25)

    return(r2_mat)
        

def calc_het(pop):
    N = float(pop.census())
    het = [C(list(pop.get_genotype(0,l).values()))[0.5]/N for l in range(pop.genomic_arch.L)]
    return(het)

    
#TODO: replace this loop with the populome approach
def calc_maf(pop):
    two_N = 2*float(pop.census())
    MAF = []
    for l in range(pop.genomic_arch.L):
        cts = C(list(pop.get_genotype(0,l).values()))
        MAF_l = (cts[1.0]*2 + cts[0.5])/two_N
        if MAF_l >0.5:
            MAF_l = 1-MAF_l
        MAF.append(MAF_l)
    return(MAF)


