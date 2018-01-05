#!/usr/bin/python

import numpy as np



#Get the phenotypic values of all individuals for a given trait
def get_phenotype(pop, trait):

    #get list of individs' ids (to make sure they are correctly associated with phenotypes
    inds = pop.individs.keys()

    #for all individuals, multiply genotype by 2*alpha (because genotype is sum of 2 alleles, each with effect size alpha) 
    #at all loci, sum across loci, then add to 0.5 (the null phenotypic value), yielding an array of all individuals' phenotypes
    phenotypes = [0.5 + sum(pop.individs[i].genome[pop.genomic_arch.traits[trait].loci,:].sum(axis = 1)*(2*pop.genomic_arch.traits[trait].alpha)) for i in inds]

    #TODO: THIS SEEMS TO BE PRODUCE STRANGELY RIGHT OR LEFT SHIFTED HISTOGRAMS OF PHENOTYPES WHEN I EXPECTED THEM
    #TO BE CENTERED ON 0.5... GO THROUGH THE MATH SLOWLY AGAIN TO SEE IF THERE'S A BUG...

    z = dict(zip(inds,phenotypes))
    return(z)



#TODO: FINISH WRITING AND TEST THIS
#NOTE: I HAVE A FEELING IF I THINK ABOUT THIS FUNCTION AND THE SERIES OF FUNCTIONS IN THIS MODULE THAT IT FITS
      #INTO, THERE WOULD BE A GOOD WAY TO SPEED IT UP
def get_fitness(pop):
    fit_dict = dict(zip(pop.individs.keys(),[1]*pop.census()))
    for t in pop.genomic_arch.traits.keys():
        s = pop.genomic_arch.traits[t].s
        gamma = pop.genomic_arch.traits[t].fitness_fn_gamma
        scape_num = pop.genomic_arch.traits[t].scape_num
        z = get_phenotype(pop, t)
        hab = pop.get_habitat(scape_num)
        w = dict([(i, 1 - s*np.abs(hab[i] - z_val)**gamma) for i,z_val in z.items()])
        fit_dict = {i : w_val *fit_dict[i] for i, w_val in w.items()}
    return(fit_dict)



#Get the vector of mortalities (probabilies of death) for a given density-dependent Pr(death) at a cell, the
#environmental value(s) at that cell, the phenotype(s) of the trait(s) for the individuals found there, and the selection
#coefficient(s) on the trait(s)
def get_prob_death(pop, d):

    W = get_fitness(pop)

    #NOTE: FOR NOW, ASSUMING d IS IN THE FORM OF A DICT OF K-V PAIRS: {indvidual_id: d_xy}, LIKE W 

    d_ind = {i: 1-(1-d[i])*w for i, w in W.items()}


    #Now, tweak it if need be to ensure that it stays between 0 and 1
    #NOTE: NEED TO CONSIDER IF THERE IS A BETTER, LESS BIASING WAY OF DOING THIS
    d_ind = {i: max(min(0.999999, d_val), 0.000001) for i, d_val in d_ind.items()}
   
    #NOTE: NEED TO FIURE OUT BEST WAY TO KEEP MORTALITIES BETWEEN 0 AND 1
    #Check that 0<= mort <= 1
    assert np.alltrue(np.array(d_ind.values())>=0)
    assert np.alltrue(np.array(d_ind.values())<=1)

    #return
    return(d_ind)

    



