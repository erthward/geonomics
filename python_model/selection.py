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
from scipy.spatial import cKDTree





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




def get_fitness(pop, norm_fit = False):  #If norm_fit == True, will normalize fitness values with respect to the max fitness value; otherwise, will leave them expressed with respect to the magnitude of the selection coefficients


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



            #filter list of each locus-respective selective environmental variable (i.e. raster numbers) associated with non-neutral loci
            env_var_non_neut = pop.genomic_arch.env_var[c][pop.genomic_arch.non_neutral[c]]
            #and use that to generate a list of each locus' current habitat value for its selective environmental variable
            hab = np.array([ind.habitat[env_var_non_neut[l] - 1] for l in range(len(env_var_non_neut))])  
            #NOTE: FOR NOW, I subtract one from env_var values to match to raster numbers for now; eventually NEED TO FIGURE OUT HOW TO MAKE 0 BE GLOBALLY SELECTIVE AND >0 BE LOCALLY SELECTIVE ACCORDING TO THE RESPECTIVE RASTER

            #filter list of selection coefficients associated with non-neutral loci
            s = np.array(pop.genomic_arch.s[c][pop.genomic_arch.non_neutral[c]])


            

            #now use the lists of non-neutral genotypes, habitat values, and selection coefficients to calculate a list of locus-specific fitnesses
            fit = s * (1 - abs(hab - g))

            
            #sum the foregoing list of locus-specific fitnesses and add it to the individual's genome-wide rolling sum of locus-specific fitnesses
            sum_fit[i] = sum_fit[i] + sum(fit)

    if norm_fit == True:
        sum_fit = normalize_fitness(pop, fit)


    assert False not in [(fitness <= 1 and fitness >= 0) for fitness in sum_fit.values()], 'Something went wrong! All individuals\' net fitness values should be between 0 and 1. Instead, returned:\n\n\t%s' % str(sum_fit)


    return sum_fit




def normalize_fitness(pop, fit):
    #divide each individual's summed locus-specific fitnesses by genomic_arch.sum_s, yielding each individual's normalized, weighted-mean fitness
    norm_fit = dict([(ind, f/pop.genomic_arch.sum_s) for ind, f in fit.items()])

    assert False not in [(fitness <= 1 and fitness >= 0) for fitness in norm_fit.values()], 'Something went wrong! All individuals\' net fitness values should be between 0 and 1. Instead, returned:\n\n\t%s' % str(norm_fit)

    return norm_fit




def select(pop, t, params, sigma_deaths = 0.05, density_dependent_fitness = True):

    #number individuals exceeding 'carrying capacity' at this timestep (as dictated by pop.size)
    overshoot = pop.census() - (pop.size[t] * pop.initial_size)

    #draw number of deaths this timestep from a normal distribution centered on the overshoot
    num_deaths = int(np.round(r.normal(overshoot, max(0.001, abs(sigma_deaths*overshoot)))))


    #get individuals' fitnesses
    fitness = get_fitness(pop)
    max_fitness = max(fitness.values())

    #express as relative negative fitnesses (ratio of difference from max to max difference from max)
    neg_fitness = [max_fitness - f for f in fitness.values()]
    #NOTE: IS THIS NEXT STEP, WHICH CONSTRAINS ALL FITNESSES TO 0 < fit < 1, JUSTIFIED, GIVEN THAT THEY ARE OFTEN DRAMATICALLY DISPARATE AT THE BEGINNING AND NARROWLY SPREAD LATER ON?...
    rel_neg_fitness = [max(0.01, min(0.99, f/max(neg_fitness))) for f in neg_fitness]






    #NOTE: TRYING TO IMPLEMENT DENSITY-DEPENDENT MORTALITY!
    if density_dependent_fitness == True:
        #create list of individs to match against list of points and list of neighbors
        individs = pop.get_coords().keys()
        #create list of points (as coords)
        points = np.array(pop.get_coords().values())
        #create cKDTree for the points
        tree = cKDTree(points, leafsize = 100)
        #create list of the count of neighbors within 5*mu_distance of each point, excluding both self and population census size (which indicates no available neighbors)
        #NOTE: The 5*mu_distance threshold is an arbitrary value for now, based on the potential combined movement of each individual in a pair 2.5 standard deviations in opposite and approaching directions...
        neighbor_counts = [len([neigh for neigh in i if (neigh != pop.census() and neigh != n)]) for n, i in enumerate(tree.query(points, k = pop.census(), distance_upper_bound = 5*params['mu_distance'])[1])]
        
        #now calculate a 'density fitness', as 1 minus the density of an individual relative to max density, and factor it directly into rel_neg_fitness (NOTE: density is never actually calculated as density proper, but this is irrelevant, as all individuals' counts would be divided by the same area (pi*(5*mu_distance)^2)  )
        #print list(rel_neg_fitness)[:10]
        rel_neg_fitness = (rel_neg_fitness + (1 - (np.array(neighbor_counts)/float(max(neighbor_counts))))) / 2.
        #print list(rel_neg_fitness)[:10]
        #raw_input()





    
    
    
    deaths = []

    while len(deaths) < num_deaths:


        dead_inds = np.array(pop.individs.keys())[np.array([bool(i) for i in r.binomial(1, rel_neg_fitness)])]

        deaths.extend(list(dead_inds))

        deaths = list(set(deaths))

    #if len(deaths) > num_deaths:

        #deaths = deaths[:num_deaths]


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


    





