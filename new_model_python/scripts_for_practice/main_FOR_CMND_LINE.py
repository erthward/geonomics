#!/usr/bin/python
#main.py


#------------------------------------------------#



#-----------------#
# IMPORT PACKAGES #
#-----------------#

import landscape
import movement
import genome
import individual
import population
import mating
import recombination
import gametogenesis
import selection
import mutation

reload(landscape)
reload(movement)
reload(genome)
reload(individual)
reload(population)
reload(mating)
reload(recombination)
reload(gametogenesis)
reload(selection)
reload(mutation)


import numpy as np
from numpy import random as r
import random
import matplotlib as mpl
import os




#------------------------------------------------#                         


#-----------------#
# GRAB PARAMETERS #
#-----------------#


#grab (and process some) command-line arguments, setting majority of "secondary" params to their default
#values ("typical", "characteristic", or uniformly-distributed values) unless otherwise stipulated in the
#command-line args






runtime = 500





#########################
      ###################
# TEMP PARAMS FOR TESTING
      ###################
#########################






params = {


'set_seed' : True,                  #set the seed (for reproducibility)?

'seed_num' : 11,                    #number to seed random number generators

'T' : runtime,                      #total model runtime

'burn_T': 100,                     #total burn-in runtime

'L' : 1e5,                         #total number of loci

'n' : 1,                           #number of chromosomes

'x' : 2,                         #ploidy (for now, leave at 2 for diploidy)

'mu' : 10e-9,                    #genome-wide mutation rate

'alpha_D' : 7e2,                #alpha for beta distribution of linkage values

'beta_D' : 7e3,                 #beta for beta distribution of linkage values


'N' : 20,                        #total pop size

'dims' : (1000,1000),             #dimensions of landscape  

'num_scapes' : 1,               #number of landscapes desired

'n_rand_pts' : 200*5,           #number of random coordinates to be used in generating the random landscapes


'interp_method' : ['linear'],   # list of interpolation methods for generation of random landscapes, 1 per landscape to be generated (as set by num_scapes)

'move' : True,                     #is this a mobile species?

'mu_direction' : 0,                #mu for von mises distribution defining movement directions

'kappa_direction' : 0,             #kappa for von mises distribution

'mu_distance' : 3,               #mean movement-distance (lognormal distribution)

'sigma_distance' : 0.1,            #sd of movement distance

'sex' : True,                      #is this a sexual species?

'repro_age' : (8, 5),          #age at sexual maturity (int or float for non-sexual species, tuple or list of two ints/floats for sexual species; set to 'None' to not make this an age-structured species


'mating_radius' : 250,              #radius of mate-searching area

'mu_dispersal' : 5,           #mean dispersal distance (lognormal distribution)

'sigma_dispersal' : 0.4,          #sd of dispersal distance

'size' : 1,              # float/int, or list/tuple of length T containing floats/ints, expressing the target population size over model time as a ratio of the starting size (N)

'sigma_deaths' : 0.2,              # std for the normal distribution used to choose the r.v. of deaths per timestep (mean of this distribution is the overshoot, as calculated from pop.size and pop.census())

'alpha_mut_s' : 25,                # alpha param for the beta distribution describing the highly advantageous selection coeffs for mutations

'beta_mut_s' : 0.5                # beta param for the beta distribution describing the highly advantageous selection coeffs for mutations



}





#########################
      ###################
      ###################
#########################


#set seed, if requested
if params['set_seed']:
    random.seed(params['seed_num'])
    r.seed(params['seed_num'])




#refresh mutation log
os.system('touch ./mutation_log.txt')



#------------------------------------------------#                         


#---------------------#
# ESTABLISH LANDSCAPE #
#---------------------#



#if landscape/s to be pulled from file, do so




#else, generate random landscapes


    
land = landscape.build_scape_stack(params = params)

  
#------------------------------------------------#                         
                                                                                                                                                               
#-----------------------------#
# CREATE GENOMIC ARCHITECTURE #
#-----------------------------#

genomic_arch = genome.build_genomic_arch(params = params, land = land)



                                                                                                                                                               
#------------------------------------------------#



#----------------------#
# ESTABLISH POPULATION #
#----------------------#



pop = population.create_population(genomic_arch = genomic_arch, land = land, params = params)





#------------------------------------------------#                         
     
     

#plot starting pop
                                                                                                                                                              
#---------#
# BURN IN #
#---------#


def burn_in(land, pop, params):

    print '\n\nSTARTING BURN-IN.\n\t(Will run for %i timesteps.)\n\n' % params['burn_T']

    pop.show(land = land, colorbar = True)

    for burn_t in range(params['burn_T']):
    
        pop.birthday()
    
        pop.move(land = land, params = params)
    
        pop.mate(land = land, params = params)

        pop.select(t = burn_t, params = params)
    
        pop.mutate(params = params, t = burn_t)
    
        pop.check_extinct()

        
        print '\n\n%i timesteps run. Current status:\n\n\t%i individuals\n\n' % (burn_t+1, pop.census())

        print '\n----------------------------------------------------------\n\n'

    
    print '\n\nBURN-IN FINISHED.\n\n'

    pop.show(land = land, colorbar = False)

    #mpl.pyplot.close()




########for t in range(1, burn_T):
########    if move:
########
########        [movement.move(ind, land, mu_movement = mu_movement, sigma_movement = sigma_movement) for ind in pop.individs.values()];
########        
########
########    
########    #find mates
########    mating_pairs = pop.find_mates(mating_radius)
########
########
########    #mate
########
########    #cross-over and produce gametes
########
########    for pair in mating_pairs:
########        #create zygote from the mating event
########        zygote = mating.mate(pop, pair, genomic_arch)
########
########
########        #generate new offspring
########
########
########        #disperse offspring
########    
########    #mutate
########    
########    #select
########
########    #add to ages
########    pop.birthday()
########
########
########
#########------------------------------------------------#                         
########                                                                                                                                                               
#########-----------#
######### RUN MODEL #
#########-----------#
########
#########once the starting, burned-in population is created, run the model for each of N different parameter sets in the input file
########
########
########for N in param_sets:
########
########    for t in range(1, T):
########        if move:
########            pass
########            # move
########
########            #find mates
########            #mate
########            #cross-over
########            #disperse offspring
########            #mutate
########            #select




#------------------------------------------------#
                                                                                                                                                               
#-----------------------------#
# POST-PROCESS RESULTING DATA #
#-----------------------------#

#post-process resulting data, to output in desired format (as specified in command-line/param-file args

