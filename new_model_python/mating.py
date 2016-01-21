#!/usr/bin/python
#mating.py

'''
##########################################

Module name:                mating

Module contains:
                            - function for simulating the mating of two individuals, according to input parameters
                            - associated functions


Author:                     Drew Ellison Hart
Email:                      drew.hart@berkeley.edu
Github:                     URL
Start date:                 12-28-15
Documentation:              URL


##########################################
'''



#TODO:

    # implement sexual selection somehow??

    # restrict males to mate with only one female per generation??





from scipy.spatial import cKDTree
import numpy as np
import numpy.random as r


import gametogenesis
import genome


#--------------------------------
# CLASSES -----------------------
#--------------------------------





#----------------------------------
# FUNCTIONS -----------------------
#----------------------------------

def find_mates(pop, mating_radius, sex = True, repro_age = None):



    #NOTE: In an IDEAL world, this would work as follows:
        #for all individuals, find set of nearest neighbors (excluding selves)
        #for females:
            #choose male from among those nn's, based on:
                #probabilistic function of distance
                #if sexual_selection:
                    #(product of?) selection coefficients




    ######################################################
    #First, build cKDTree and create nearest-neighbor query:

    individs = [(i, [ind.x, ind.y], ind.sex, ind.age) for i, ind in pop.individs.items()]
    points = np.array([ind[1] for ind in individs])
    tree = cKDTree(points, leafsize = 100)   #NOTE: Figure out how leafsize, and how to parameterize this in order to optimize speed for any given pop...
    query = tree.query(points, k = 2, distance_upper_bound = mating_radius)
    #query_2 = tree.query_ball_tree(other = tree, r = mating_radius) #alternative, using query_ball_tree() method
    #[query_2[i].remove(i) for i in range(len(query_2))]; #remove self-pairings




    ####################################################
    # Then, operationalize sexes, if being used, and find all available pairs within max distance

    if sex:
        #np.array of the sexes of all individuals
        sexes = np.array([ind[2] for ind in individs])


        #array of couplings for all females with nearest individual < mating_radius
        available_females = np.array(query[0][:,1] != np.inf) & np.array(sexes[query[1][:,0]] == 0) # i.e.  [AT LEAST 1 INDIVID < mating_radius (OTHERWISE scipy.spatial.cKDTree associates them with a nonexistent neighbour and infinite distance)] & [female]

        #check that nearest individuals to paired females are males, otherwise remove
        available_pairs = query[1][available_females][sexes[query[1][:,1][available_females]] == 1]
        available_pairs_inds = [query[1][i] in available_pairs for i in range(len(query[1]))]


    else:

        #take all pairs within mating_radius
        available_individs = np.array(query[0][:,1] != np.inf) 

        #delete inverse-equal pairs from list (i.e. if 3-paired-to-5 and 5-paired-to-3 are both listed, drop one)
        available_pairs = np.array([np.array(tup) for tup in list(set([tuple(set(item)) for item in query[1][available_individs]]))])
        available_pairs_inds = [query[1][i] in available_pairs for i in range(len(query[1]))]






    ###############################################
    # Then, operationalize age-structuing, if being used


        #NOTE: I have NOT restricted this so that each individual can only mate once per generation, but I should probably figure that out and add it here!

    if repro_age <> None:
        #np.array of the ages of all individuals
        ages = np.array([ind[3] for ind in individs])

        if sex: #if sexual species, repro_age expected to be a tuple or list of numerics of length 2

            assert repro_age.__class__.__name__ in ('ndarray', 'list', 'tuple'), "For a sexual and age-structured species, age at sexual maturity, 'repro_age', must be expressed as an iterable of numerics (i.e. a list, tuple or numpy.ndarray of floats or integers); you have provided a:\n\t%s" % repro_age.__class__.__name__


            yes_f = np.array(ages[available_pairs[:,0]] >= repro_age[0])
            yes_m = np.array(ages[available_pairs[:,1]] >= repro_age[1])
            available_pairs = available_pairs[yes_f & yes_m]
            available_pairs_inds = [query[1][i] in available_pairs for i in range(len(query[1]))]

        else: #if not sexual species, repro_age expected to be a numeric

            assert repro_age.__class__.__name__ in ('int', 'float'), "For a non-sexual and age-structured species, the age at sexual maurity, 'repro_age', must be expressed as a single, non-iterable numeric (i.e.  float or integer); you have provided a:\n\t%s" % repro_age.__class__.__name__


            yes_0 = np.array(ages[available_pairs[:,0]] >= repro_age)
            yes_1 = np.array(ages[available_pairs[:,1]] >= repro_age)
            available_pairs = available_pairs[yes_0 & yes_1]
            available_pairs_inds = [query[1][i] in available_pairs for i in range(len(query[1]))]




    #########################################
    # Then, if at least one valid pair was found, decide whether or not it will reproduce, as a binomial random var weighted by pair distance(s)


    #if valid pairs exist:
    if len(available_pairs) > 0:

        #binomial decision whether or not to mate with nearest male, as function of ratio: nn distance/mating_radius 
        mating_decisions = np.array([r.binomial(n = 1, p = p) for p in 1-(query[0][available_pairs[:,0]][:,1]/mating_radius) ]) == 1

        mates = available_pairs[mating_decisions]


        #finally, link back to initially created structure, to get individuals' proper keys
        keys = [i[0] for i in individs]

        mates = np.array([[keys[mate] for mate in pair] for pair in mates])
        

    
    else:

        mates = np.array([])

    ####################################################
    #Return an array or arrays, each inner array containing a mating-pair
    return mates








#function for mating a chosen mating-pair
def mate(pop, pair, genomic_arch, n_num_offspring = 1, p_num_offspring = 0.7):
    offspring = []
    num_offspring = r.negative_binomial(n_num_offspring, p_num_offspring) + 1
    for n in range(num_offspring):
        zygote = {}
        gametes = [gametogenesis.gametogenerate(pop.individs[i], genomic_arch) for i in pair]
        for c in range(genomic_arch.n):
            chromosome = np.vstack((gametes[0][c], gametes[1][c])).T
            zygote[c] = chromosome

        offspring.append(genome.Genome(zygote, genomic_arch.x))

    return offspring




