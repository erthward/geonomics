#!/usr/bin/python
#selection.py


'''
##########################################

Module name:                selection

Module contains:
                            - function for simulating selection at a single locus, according to input parameters
                            - functions for simulating epistatic and/or polygenic selection
                            - associated functions


Author:                     Drew Ellison Hart
Email:                      drew.hart@berkeley.edu
Github:                     URL
Start date:                 12-28-15
Documentation:              URL


##########################################
'''




##########

#TODO:
    #some means of allowing selection at different loci for different rasters

    #additive selection? genic selection?

    #dominance? codominance?

    #loci that are universally selective? loci that are universally selective, but most strongly in a 1 or 0 habitat (i.e raster cell)? loci whose two alleles are resectively adapative to 0 and 1 habitat (i.e.raster cell values)?

    #directional selection? hetero advantage? hetero disadvantage?

    #fertility selection AND/OR viability selection? 

    #sexual selection?
                   

##########




import numpy as np
from numpy import random as r





#------------------------------------
# CLASSES ---------------------------
#------------------------------------



#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------

#
#alpha = 3
#
#beta = 0.3
#
#s_dist = lambda: r.beta(alpha, beta)
#
#coeff = s_dist()
#
#locus = np.array(range(len(genomic_arch.s[0])))[genomic_arch.s[0] == min(genomic_arch.s[0], key = lambda x: abs(x-coeff))][0]
#
#habitats = pop.return_habitat()
#
#genotypes = pop.return_genotypes(locus, chromosome_num = 0, format = 'mean')
#
#selection_strengths = [(i, 1 - abs(habitats[i] - genotypes[i])) for i in habitats.keys()]




def get_fitness(pop):


    #dict of functions correpsonding to the various dominance types
    dom_funcs = { 0 : np.mean,                          #codominance
                  1 : lambda x: np.ceil(np.mean(x)),    #dominant (with respect to relationship between allele 1 and habitat values -> 1
                  2 : lambda x: np.floor(np.mean(x))    #recessive (with respect to relationship between allele 1 and habitat values -> 1
                }


    #initiate dict that will eventually hold all individuals' net fitnesses
    sum_fit = dict(zip(pop.individs.keys(), [0]*pop.census()))


    #for each chromosome (seems strange, but programatically cleaner to make this the outer loop)
    for c in pop.genomic_arch.s.keys():


        #for each individual
        for i, ind in pop.individs.items():


            #filter list of dominance types associated with non-neutral loci
            dom_non_neut = pop.genomic_arch.dom[c][pop.genomic_arch.non_neutral[c]]
            #filter list of genotypes associated with non-neutral loci
            g_non_neut = ind.genome.genome[c][pop.genomic_arch.non_neutral[c]]
            #and use those two lists to calculate the genotypes for selective loci in 0.0/0.5/1.0-format
            g = np.array([dom_funcs[dom_non_neut[l]](g_non_neut[l]) for l in range(len(dom_non_neut))])



            #filter list of each locus-respective selective environmental variables (i.e. raster numbers) associated with non-neutral loci
            env_var_non_neut = pop.genomic_arch.env_var[c][pop.genomic_arch.non_neutral[c]]
            #and use that to generate a list of each locus' current habitat value for its selective environmental variable
            hab = np.array([ind.habitat[env_var_non_neut[l] - 1] for l in range(len(env_var_non_neut))])  
            #NOTE: FOR NOW, I subtract one from env_var values to match to raster numbers for now; eventually NEED TO FIGURE OUT HOW TO MAKE 0 BE GLOBALLY SELECTIVE AND >0 BE LOCALLY SELECTIVE ACCORDING TO THE RESPECTIVE RASTER

            #filter list of selection coefficients associated with non-neutral loci
            s = np.array(pop.genomic_arch.s[c][pop.genomic_arch.non_neutral[c]])


            

            #now use the lists of non-neutral genotypes, habitat values, and selection coefficients to calculate a list of locus-specific fitnesses
            fit = s * (1 - abs(hab - g))

            
            #sum the foretold list of locus-specific fitnesses and add it to the individual's genome-wide rolling sum of locus-specific fitnesses
            sum_fit[i] = sum_fit[i] + sum(fit)



    #finally, divide each individual's summed locus-specific fitnesses by genomic_arch.sum_s, yielding each individual's normalized, weighted-mean fitness
    sum_fit = dict([(ind, fit/pop.genomic_arch.sum_s) for ind, fit in sum_fit.items()])

    assert False not in [(fitness <= 1 and fitness >= 0) for fitness in sum_fit.values()], 'Something went wrong! All individuals\' net fitness values should be between 0 and 1. Instead, returned:\n\n\t%s' % str(sum_fit)

    return sum_fit




def select(pop, t, sigma_deaths = 0.2):

    #number individuals exceeding 'carrying capacity' at this timestep (as dictated by pop.size)
    overshoot = pop.census() - (pop.size[t] * pop.initial_size)

    #draw number of deaths this timestep from a normal distribution centered on the overshoot
    num_deaths = int(np.round(r.normal(overshoot, sigma_deaths)))



    #get individuals' fitnesses
    fitness = get_fitness(pop)

    #express as relative negative fitnesses (ratio of difference from max to max difference from max)
    neg_fitness = [max(fitness.values()) - f for f in fitness.values()]
    rel_neg_fitness = [max(0.01, min(0.99, f/max(neg_fitness))) for f in neg_fitness]

    deaths = []

    while len(deaths) < num_deaths:


        dead_inds = np.array(pop.individs.keys())[np.array([bool(i) for i in r.binomial(1, rel_neg_fitness)])]

        deaths.extend(list(dead_inds))

        deaths = list(set(deaths))

    if len(deaths) > num_deaths:

        deaths = deaths[:num_deaths]


    print '\t%i individuals dead' % len(deaths)
    
    [pop.individs.pop(ind) for ind in deaths]






####NOTES, PREVIOUS

    #universal selection


    
    #spatially variable selection, heterozygote advantage in intermediate habitat (i.e. codominance?)



    #spatially variable selection, dominant homozygote advantage in dominant's habitat (i.e. dominance?)






    #NOTE: OR SHOULD I JUST KILL OFF HOWEVER MANY DIE BASED ON FITNESS, AND THEN LET THE POP SIZE CONTROL BE IMPLEMENTED DURING THE BREEDING PERIOD???



    #compose some function for compiling an individual's overall fitness based on its genotype at various loci



    #vectorize a calculation based on these composite fitnesses, and make those whose composite fitnesses match their habitat values more closely more fit



    #have some way to operationalize UNIVERSALLY advantageous alleles separately from LOCALLY advantageous alleles!!


    





