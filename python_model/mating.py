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

#import sklearn.neighbors  #added by Irene

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

def find_mates(pop, params, sex = False, repro_age = None, dist_weighted_birth = False):

    b = params['b']
    mating_radius = params['mating_radius']
    if 'sex' in params.keys():
        sex = params['sex']
    if 'repro_age' in params.keys():
        repro_age = params['repro_age']
    



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
    ##query_2 = tree.query_ball_tree(other = tree, r = mating_radius) #alternative, using query_ball_tree() method
    ##[query_2[i].remove(i) for i in range(len(query_2))]; #remove self-pairings

    #nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(points) #added by irene
    #distances, indices = nbrs.kneighbors(X) #added by irene




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




    ##############################################
    #Check if any available pairs have been found thus far, proceed if so, otherwise return empty array of pairs
    if len(available_pairs) == 0:
        mates = np.array([])
        return mates
    else:
        pass



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

        #if birth probability is not to be weighted by the distance between individuals in a pair
        if dist_weighted_birth == False:

            #binomial decision whether or not to mate with nearest male, as function of ratio: nn distance/mating_radius 
            mating_decisions = r.binomial(n = 1, p = b, size = available_pairs.shape[0]) == 1

            mates = available_pairs[mating_decisions]


            #finally, link back to initially created structure, to get individuals' proper keys
            keys = [i[0] for i in individs]

            mates = np.array([[keys[mate] for mate in pair] for pair in mates])

        #if birth probability is to be weighted by the distance btwn individs in a pair
        elif dist_weighted_birth == True:
 
            #binomial decision whether or not to mate with nearest male, as function of ratio: nn distance/mating_radius 
            mating_decisions = np.array([r.binomial(n = 1, p = b) for p in 1-(query[0][available_pairs[:,0]][:,1]/mating_radius) ]) == 1

            mates = available_pairs[mating_decisions]


            #finally, link back to initially created structure, to get individuals' proper keys
            keys = [i[0] for i in individs]

            mates = np.array([[keys[mate] for mate in pair] for pair in mates])

       

    
    else:

        mates = np.array([])

    ####################################################
    #Return an array or arrays, each inner array containing a mating-pair
    return mates






def determine_num_births(num_pairs, params):
    lambda_offspring = params['lambda_offspring']
    num_births = [r.poisson(lambda_offspring-1+.0001) + 1 for i in range(num_pairs)]
    return(num_births)



#function for mating a chosen mating-pair
def mate(pop, pair, params, n_offspring, gamete_recomb_paths):
    offspring = []
        #NOTE: Subtracting 1 from the input lambda then
        #adding 1 to the resulting vector ensures that all pairs who were already determined to be giving birth
        #will have at least one offspring (i.e. prevents draws of 0); adding .000001 to the input lambda prevents
        #absolute fixation at 1 that results when input lambda is 0
        #NOTE: This means that the deaths are actually not truly Poisson dispersed (actually underdispersed, I
        #believe); I don't really see that this matters much, given that I don't even want to keep this simply
        #a Poisson dist in the long-run (would rather change to neg binom, to allow overdispersion, with option
        #to parameterize it to be equivalent to a Poisson), but even still I would need to use the
        #(lambda-1) + 1 trick either way to avoid 0s, so that issue would remain; for now, minor, but should
        #come back and think the theoretical implications/justification for this
    for n in range(n_offspring):
        #generate a gamete for each member of mating pair
        #paths, gamete_recomb_paths = [gamete_recomb_paths[:,0], gamete_recomb_paths[:,1]], gamete_recomb_paths[:,2:]
        paths = [gamete_recomb_paths.pop() for i in range(2)]
        gametes = [gametogenesis.gametogenerate(pop.individs[i], paths.pop()) for i in pair]
        #stack the gametes and transpose, to create the new individual's new_genome array
        new_genome = np.vstack((gametes[0], gametes[1])).T
        #append the new_genome to the list of offspring
        offspring.append(new_genome)

    return offspring




