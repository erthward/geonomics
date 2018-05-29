#!/usr/bin/python

import numpy as np
from collections import OrderedDict as OD
from itertools import repeat


#Get the phenotypic values of all individuals for a given trait
def calc_phenotype(ind, genomic_arch, trait):

    #get number of loci and their allelic effect sizes
    n_loci = trait.n_loci
    alpha = trait.alpha

    #if monogenic, then divide genotype by 2, 
    #but if polygenic then multiply genotype by alpha (i.e. effect size) 
    #at all loci, sum across loci, then add to 0.5 (the null phenotypic value)
    phenotype = np.float32(ind.genome[trait.loci,:].sum(axis = 1))
    if n_loci > 1:
        phenotype = 0.5 + sum(phenotype*alpha)
    else:
        phenotype = phenotype[0]/2

    return(phenotype)



def get_single_trait_fitness(pop, trait):
    trt = pop.genomic_arch.traits[trait]
    scape_num = trait.scape_num
    hab = pop.get_habitat_by_land_ind(scape_num = scape_num)
    phi = trait.get_phi(pop)
    gamma = trait.gamma
    scape_num = trait.scape_num
    z = pop.get_phenotype()
    univ_advant = trt.univ_advant
    w = {i:1 - phi[i]*(abs((hab[i]**(not univ_advant)) - z_val[trait])**gamma) for i,z_val in z.items()}
    return(w)

def calc_sngl_trt_fitness(t, e, z, pop):
    fit = 1 - t.get_phi(pop)*(abs((e[:,t.scape_num]**(not t.univ_advant)) - z[:,t.num])**t.gamma)
    return(fit)

def calc_fitness(pop, traits):
    #get all individuals' environmental values
    e = pop.get_habitat()
    #get all individuals' phenotypes
    z = pop.get_phenotype()
    #create lambda function with current e, z, and pop objects
    calc_fitness_lambda = lambda t: calc_sngl_trt_fitness(t, e, z, pop)
    #map the calc_sngl_trait_fitness function to all traits
    fit = np.stack(list(map(calc_fitness_lambda, traits))).prod(axis = 0)
    return(fit)


#Get the vector of mortalities (probabilies of death) for a given density-dependent Pr(death) at a cell, the
#environmental value(s) at that cell, the phenotype(s) of the trait(s) for the individuals found there, and the selection
#coefficient(s) on the trait(s)
def get_prob_death(pop, d):

    W = get_fitness(pop)

    #NOTE: FOR NOW, ASSUMING d IS IN THE FORM OF A DICT OF K-V PAIRS: {indvidual_id: d_xy}, LIKE W 

    death_probs = OD()

    [death_probs.update({i: 1-(1-d[i])*w}) for i, w in W.items()];
   
    assert np.alltrue(np.array(list(death_probs.values()))>=0)
    assert np.alltrue(np.array(list(death_probs.values()))<=1)

    assert False not in [w >= 0 and w <= 1 for w in death_probs.values()], 'ERROR: Some fitness values outside the 0-to-1 range.'

    #return
    return(death_probs)

    



