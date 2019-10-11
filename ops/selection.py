#!/usr/bin/python
#selection.py


'''
##########################################

Module name:          ops.selection


Module contains:
                      - Functions to implement fitness and selection


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           04-28-17
Documentation:        URL


##########################################
'''


import numpy as np
from collections import OrderedDict as OD


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#Get the phenotypic values of all individuals for a given trait
def _calc_phenotype(ind, gen_arch, trait):
    #get number of loci and their allelic effect sizes
    n_loci = trait.n_loci
    alpha = trait.alpha
    #get the mean genotype array
    genotype = np.mean(ind.g[trait.loci], axis = 1)
    #use dominance, if required (to save considerable compute time otherwise)
    if gen_arch._use_dom:
        #get the dominance values
        dom = gen_arch.dom[trait.loci]
        #update the genotype array by accounting for dominance at each locus
        genotype = np.clip(genotype * (1 + dom), a_min = None, a_max = 1)
    #if polygenic, multiply genotype by alpha (i.e. effect size) 
    #at all loci, sum across loci, then add to 0.5 (the null phenotypic
    #value), to get phenotype
    if n_loci > 1:
        phenotype = 0.5 + sum(genotype*alpha)
    #else if monogenic, then mean genotype = phenotype
    else:
        phenotype = genotype[0]
    return(phenotype)


def _calc_fitness_one_trait(t, e, z, spp):
    fit = 1 - t._get_phi(spp)*(abs((
                e[:,t.lyr_num]**(not t.univ_adv)) - z[:,t.idx])**t.gamma)
    return(fit)


def _calc_fitness_traits(spp, trait_num = None):
    traits = spp.gen_arch.traits.values()
    #subset for single trait, if indicated
    if trait_num is not None:
        traits = [list(traits)[trait_num]]
    #get all individuals' environmental values
    e = spp._get_e()
    #get all individuals' phenotypes
    z = spp._get_z()
    #create lambda function with current e, z, and spp objects
    calc_fitness_lambda = lambda t: _calc_fitness_one_trait(t, e, z, spp)
    #map the calc_sngl_trait_fitness function to all traits, then
    #calculate overall fitness as product of
    #fitness for each trait
    fit = np.stack(list(map(calc_fitness_lambda, traits))).prod(axis = 0)
    #for polygenic traits, loci with very large effect sizes can 
    #generate fitnesses less than 0; correct for this
    fit = np.clip(fit, a_min = 0.001, a_max = None)
    return(fit)


def _calc_fitness_deleterious_mutations(spp):
    #create an np.array that has individuals in rows and their
    #diploid genotypes for each of the deleterious loci in the cols
    #(0, 1, or 2, to facilitate the fitness math, because s values
    #(i.e. selection coefficients) are expressed per allele)
    deletome = np.sum(np.stack([ind.g[[*spp.gen_arch.delet_loci.keys(
                                    )],:] for ind in spp.values()]), axis = 2)
    fit = 1 - np.multiply(deletome, np.array(
                                        [*spp.gen_arch.delet_loci.values()]))
    fit = fit.prod(axis = 1)
    return(fit)


#one function to calculate total fitness, including traits and deleterious
#loci, as applicable
def _calc_fitness(spp, trait_num=None):
    #set a default w array
    w = np.array([1]*len(spp))
    #get trait-related fitness, if traits
    if (spp.gen_arch.traits is not None
        and len(spp.gen_arch.traits) > 0):
        w = w * _calc_fitness_traits(spp, trait_num = trait_num)
    #if fitness is not supposed to be calculated for a specific trait, and if 
    #species has deleterious mutations, then get the fitnesses related to 
    #the deleterious traits (i.e. operationalize background selection)
    if (trait_num is None and
        len(spp.gen_arch.delet_loci) > 0):
            w = w * _calc_fitness_deleterious_mutations(spp)
    return w


#Get the vector of mortalities (probabilies of death) for a 
#given density-dependent Pr(death) at a cell, the environmental value(s)
#at that cell, the phenotype(s) of the trait(s) for the individuals found
#there, and the selection coefficient(s) on the trait(s)
def _calc_prob_death(spp, d):
    #get the fitness values (while also setting all individ.fit attributes)
    w = spp._calc_fitness()
    death_probs = 1-(1-d)*w
    assert (death_probs >= 0).all() and (death_probs <= 1).all(), ("Some "
                        "death-probability values outside the 0-to-1 range.")
    return(death_probs)

