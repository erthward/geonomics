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

import numpy as np
from numpy import random as r
import random
import matplotlib as mpl
import os



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




#------------------------------------------------#                         


#-----------------#
# GRAB PARAMETERS #
#-----------------#


#grab (and process some) command-line arguments, setting majority of "secondary" params to their default
#values ("typical", "characteristic", or uniformly-distributed values) unless otherwise stipulated in the
#command-line args






runtime = 5000





#######################
      #################
      # SET TEMP PARAMS
      #################
#######################






params = {


'set_seed' : True,                  #set the seed (for reproducibility)?

'seed_num' : 1,                    #number to seed random number generators

'T' : runtime,                      #total model runtime

'burn_T': 50,                     #total burn-in runtime

'L' : 5e2,                         #total number of loci

'n' : 1,                           #number of chromosomes

'x' : 2,                         #ploidy (for now, leave at 2 for diploidy)

'n_traits' : 150,                 #number of traits to simulate

'mu' : 10e-9,                    #genome-wide mutation rate

'alpha_D' : 7e2,                #alpha for beta distribution of linkage values

'beta_D' : 7e3,                 #beta for beta distribution of linkage values

'use_dom' : False,              #whether or not to use dominance (default to False)
                                #NOTE: REALLY JUST NEED TO GET RID OF THE DOMINANCE THING; IT'S ALL MESSED UP

'N' : 100,                        #total pop size

'dims' : (25,25),             #dimensions of landscape  

'num_scapes' : 1,               #number of landscapes desired

'n_rand_pts' : 200*5,           #number of random coordinates to be used in generating the random landscapes

#'interp_method' : ['nearest'],
'interp_method' : ['linear'],   # list of interpolation methods for generation of random landscapes, 1 per landscape to be generated (as set by num_scapes)

'move' : True,                     #is this a mobile species?

'movement_surf' : True,       #use a landscape layer as a resistance surface, or habitat quality layer, to direct movement?
#'movement_surf' : False,

'movement_surf_scape_num' : 0,               #scape number to use as the movement surface

'movement_surf_vonmises_kappa' : 2, #kappa value to use in the von Mises mixture distributions (KDEs) underlying resistance surface movement

'movement_surf_gauss_KDE_bandwidth' : 0.2, #bandwidth value to use in the Gaussian KDEs that are created to approximate the von Mises mixture distributions (KDEs) underlying resistance surface movement

'mu_direction' : 0,                #mu for von mises distribution defining movement directions

'kappa_direction' : 0,             #kappa for von mises distribution

'mu_distance' : 0.01,               #mean movement-distance (lognormal distribution)

'sigma_distance' : 0.001,            #sd of movement distance

'sex' : False,                      #is this a sexual species?

'repro_age' : 5,          #age at sexual maturity (int or float for non-sexual species, tuple or list of two ints/floats for sexual species; set to 'None' to not make this an age-structured species


'mating_radius' : 2,              #radius of mate-searching area

'mu_dispersal' : 0.1,           #mean dispersal distance (lognormal distribution)

'sigma_dispersal' : 0.02,          #sd of dispersal distance

'size' : 1,              # float/int, or list/tuple of length T containing floats/ints, expressing the target population size over model time as a ratio of the starting size (N)

'sigma_deaths' : 0.2,              # std for the normal distribution used to choose the r.v. of deaths per timestep (mean of this distribution is the overshoot, as calculated from pop.size and pop.census())

'density_dependent_fitness' : True, #should fitness be density dependent? (note: helps to avoid subpopulation 'clumping')

'alpha_mut_s' : 25,                # alpha param for the beta distribution describing the highly advantageous selection coeffs for mutations

'beta_mut_s' : 0.5                # beta param for the beta distribution describing the highly advantageous selection coeffs for mutations



}





###############################
      #########################
      # CREATE MODEL COMPONENTS
      #########################
###############################


#set seed, if requested
if params['set_seed']:
    random.seed(params['seed_num'])
    r.seed(params['seed_num'])




#refresh mutation log
os.system('touch /home/ihavehands/Desktop/mutation_log.txt')



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
     
     

     
######################
      ################
      # RUN MAIN MODEL
      ################
######################                                                                                                                                                         


#plot starting pop

#---------#
# BURN IN #
#---------#


def burn_in(pop, land, params):

    print '\n\nSTARTING BURN-IN.\n\t(Will run for %i timesteps.)\n\n' % params['burn_T']

    #pop.show(land = land, colorbar = True)

    for burn_t in range(params['burn_T']):

        print(burn_t)

        pop.birthday()
    
        pop.move(land = land, params = params)
    
        pop.mate(land = land, params = params)

        pop.select(t = burn_t, params = params)
    
        pop.mutate(params = params, t = burn_t)
    
        pop.check_extinct()

        
        print '\n\n%i timesteps run. Current status:\n\n\t%i individuals\n\n' % (burn_t+1, pop.census())

        print '\n----------------------------------------------------------\n\n'

    
    print '\n\nBURN-IN FINISHED.\n\n'

    #pop.show(land = land, colorbar = False)

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

