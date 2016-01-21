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

reload(landscape)
reload(movement)
reload(genome)
reload(individual)
reload(population)
reload(mating)
reload(recombination)
reload(gametogenesis)
reload(selection)


import numpy as np
from numpy import random as r
import random
import matplotlib as mpl




#------------------------------------------------#                         


#-----------------#
# GRAB PARAMETERS #
#-----------------#


#grab (and process some) command-line arguments, setting majority of "secondary" params to their default
#values ("typical", "characteristic", or uniformly-distributed values) unless otherwise stipulated in the
#command-line args






#########################
      ###################
# TEMP PARAMS FOR TESTING
      ###################
#########################

set_seed = True

seed_num = 11

T = 100                         #total model runtime

L = 1e5                         #total number of loci

n = 1                           #number of chromosomes

N = 10                         #total pop size

dims = (200,200)              #dimensions of landscape  

move = True                     #is this a mobile species?

mu_direction = 0                #mu for von mises distribution defining movement directions

kappa_direction = 0             #kappa for von mises distribution

mu_distance = 3               #mean movement-distance (lognormal distribution)

sigma_distance = 0.1            #sd of movement distance

sex = True                      #is this a sexual species?

repro_age = (8, 5)           #age at sexual maturity (int or float for non-sexual species, tuple or list of two ints/floats for sexual species; set to 'None' to not make this an age-structured species


mating_radius = 250              #radius of mate-searching area

mu_dispersal = 5                #mean dispersal distance (lognormal distribution)

sigma_dispersal = 0.4          #sd of dispersal distance

size = [1] * T                  # float/int, or list/tuple of length T containing floats/ints, expressing the target population size over model time as a ratio of the starting size (N)


mu = 10e-9                     #genome-wide mutation rate


#########################
      ###################
      ###################
#########################


#set seed, if requested
if set_seed:
    random.seed(seed_num)
    r.seed(seed_num)




#------------------------------------------------#                         


#---------------------#
# ESTABLISH LANDSCAPE #
#---------------------#



#if landscape/s to be pulled from file, do so




#else, generate random landscapes


num_scapes = 2
    
land = landscape.build_scape_stack(num_scapes, dims, round(1.3*dims[0]), interp_method = ['linear', 'nearest'])

  
#------------------------------------------------#                         
                                                                                                                                                               
#-----------------------------#
# CREATE GENOMIC ARCHITECTURE #
#-----------------------------#

genomic_arch = genome.build_genomic_arch(L, n, mu, land)



                                                                                                                                                               
#------------------------------------------------#



#----------------------#
# ESTABLISH POPULATION #
#----------------------#



pop = population.create_population(N = N, genomic_arch = genomic_arch, dims = dims, size = size, land = land, T = T)




#------------------------------------------------#                         
                                                                                                                                                               
#---------#
# BURN IN #
#---------#


for t in range(10):
    pop.move(land)
    pop.mate(land, mating_radius, mu_dispersal, sigma_dispersal, sex = True)
    selection.select(pop, t)

    pop.show(land, 0)

    mpl.pyplot.close()


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

