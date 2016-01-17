#!/usr/bin/python
#main.py


#------------------------------------------------#



#-----------------#
# IMPORT PACKAGES #
#-----------------#


import genome
import population
import recombination
import gametogenesis

import numpy as np
from numpy import random as r
import random




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

seed_num = 3

sex = True                      #is this a sexual species?

repro_age = (8, 5)           #age at sexual maturity (int or float for non-sexual species, tuple or list of two ints/floats for sexual species; set to 'None' to not make this an age-structured species

move = True                     #is this a mobile species?

T = 100                         #total model runtime

L = 1e5                         #total number of loci

n = 1                           #number of chromosomes

N = 10                         #total pop size

dims = (1000,1000)              #dimensions of landscape  

mu_movement = 3               #mean movement-distance (lognormal distribution)

sigma_movement = 0.1            #sd of movement distance

mating_radius = 250              #radius of mate-searching area

mu_dispersal = 5                #mean dispersal distance (lognormal distribution)

sigma_dispersal = 0.4          #sd of dispersal distance


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




import landscape
    
    #NOTE: FOR NOW, JUST USING A SINGLE LANDSCAPE LAYER
    #landscapes = dict()
    
    #for i in range(num_landscapes):
      #landscapes[i] = gen_land(params) 

    
land = landscape.random_surface(dims,round(0.7*dims[0]))

  
#------------------------------------------------#                         
                                                                                                                                                               
#-----------------------------#
# CREATE GENOMIC ARCHITECTURE #
#-----------------------------#

genomic_arch = genome.build_genomic_arch(L, n)



                                                                                                                                                               
#------------------------------------------------#



#----------------------#
# ESTABLISH POPULATION #
#----------------------#



pop = population.create_population(N, genomic_arch, dims, 0.2, 0.2, T) 




#------------------------------------------------#                         
                                                                                                                                                               
#---------#
# BURN IN #
#---------#


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

