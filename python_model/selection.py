#!/usr/bin/python

import numpy as np



#Get the phenotypic values of all individuals for a given trait
def get_phenotype(pop, trait):

    #get list of individs' ids (to make sure they are correctly associated with phenotypes
    inds = pop.individs.keys()

    #get the trait object
    trait_obj = pop.genomic_arch.traits[trait]

    #for all individuals, multiply genotype by 2*alpha (because genotype is sum of 2 alleles, each with effect size alpha) 
    #at all loci, sum across loci, then add to 0.5 (the null phenotypic value), yielding an array of all individuals' phenotypes
        #NOTE: NO! That was a mistake to multiply alpha by 2 manually. 
        #Because when I sum across axis 1, then individuals who homozygous 1 sum
        #to 2, which is then multiplied by alpha accordingly. I think I was mistakenly 
        #thinking that I was summing mean environmental values, i.e. 0.5's
    phenotypes = [(trait_obj.n > 1)*0.5 + sum(pop.individs[i].genome[trait_obj.loci,:].sum(axis = 1)*(trait_obj.alpha)) for i in inds]
        #NOTE: DEH 01-06-17: Replaced 0.5 with (trait_obj.n >1)*0.5 in the line above, so that for multigenic
        #traits, 0.5 is treated as the mean phenotypic value to which positive or negative alpha deviations,
        #whereas for monogenic traits a positive alpha effect is added to the default phenotypic value of 0

    #TODO: THIS SEEMS TO BE PRODUCE STRANGELY RIGHT OR LEFT SHIFTED HISTOGRAMS OF PHENOTYPES WHEN I EXPECTED THEM
    #TO BE CENTERED ON 0.5... GO THROUGH THE MATH SLOWLY AGAIN TO SEE IF THERE'S A BUG...
        #NOTE: 01-06-17: NOPE, JUST BECAUSE OF SAMPLING OF ALPHAS, WHERE A FEW OUTLIERS HAVE AN OUTSIZE EFFECT ON THE PHENOTYPIC DIST

    z = dict(zip(inds,phenotypes))
    return(z)



#NOTE: I HAVE A FEELING IF I THINK ABOUT THIS FUNCTION AND THE SERIES OF FUNCTIONS IN THIS MODULE THAT IT FITS
      #INTO, THERE WOULD BE A GOOD WAY TO SPEED IT UP

#TODO: DEH 01-11-18: Added ability to use spatially varying phenotypic selection coefficient, but curious if 
#it would be faster to add a conditional check of whether or not a given trait's phi is a single value, and to
#only calculate as phi[i]*(np.abs(hab[i][scape_num]... if it isn't??
def get_fitness(pop):
    fit_dict = dict(zip(pop.individs.keys(),[1]*pop.census()))
    hab = pop.get_habitat()
    for t,trt in pop.genomic_arch.traits.items():
        phi = trt.get_phi(pop)
        gamma = trt.fitness_fn_gamma
        scape_num = trt.scape_num
        z = get_phenotype(pop, t)
        univ_advant = trt.univ_advant
        w = {i:1 - phi[i]*(np.abs((hab[i][scape_num]**(not univ_advant)) - z_val)**gamma) for i,z_val in z.items()}
        fit_dict = {i : w_val *fit_dict[i] for i, w_val in w.items()}
    return(fit_dict)



#NOTE: I HAVE A FEELING IF I THINK ABOUT THIS FUNCTION AND THE SERIES OF FUNCTIONS IN THIS MODULE THAT IT FITS
      #INTO, THERE WOULD BE A GOOD WAY TO SPEED IT UP
def get_single_trait_fitness(pop, trait):
    hab = pop.get_habitat_by_land_ind(scape_num = scape_num)
    trt = pop.genomic_arch.traits[trait]
    phi = trt.get_phi(pop)
    gamma = trt.fitness_fn_gamma
    scape_num = trt.scape_num
    z = get_phenotype(pop, trait)
    univ_advant = trt.univ_advant
    w = {i:1 - phi[i]*(np.abs((hab[i]**(not univ_advant)) - z_val)**gamma) for i,z_val in z.items()}
    return(w)



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

    



