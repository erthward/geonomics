#!/usr/bin/python
#main.py


'''
    @@@@@  @@@@@  @@@@@       @    @  @@@@@  @     @  @@@@@  @@@@@  @@@@@
   @      @      @   @       @@   @  @   @  @@   @@    @    @      @    
  @  @@  @@@@@  @   @  @@@@ @ @  @  @   @  @ @ @ @    @    @      @@@@@
 @   @  @      @   @       @  @ @  @   @  @  @  @    @    @          @
@@@@@  @@@@@  @@@@@       @    @  @@@@@  @     @  @@@@@  @@@@@  @@@@@


GEO-NOMICS: A generalizable, spatially explicit, 
            individual-based, forward-time 
            landscape-genomic simulation model



Author:         Drew Ellison Hart
Email:          drew.hart@berkeley.edu
Github:         URL
Start date:     12-28-15
Documentation:  URL



Usage:

  IMPORTANT NOTE:
      docs should include an example of a typical param-file, and then explain that additional params (i.e.
      params not expected to be commonly tweaked/explored/analyzed) could be
      included therein, or as command-line args, in order to alter them from their default settings
  

'''



#if -h or --help flag in command-line args, return docstring and exit
if '-h' in sys.argv or '--help' in sys.argv:
    print(__doc__)
    sys.exit()



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

T = 100                         #total model runtime

L = 1e5                         #total number of loci

n = 1                           #number of chromosomes

N = 100                         #total pop size

dims = (1000,1000)              #dimensions of landscape  

move = True                     #is this a mobile species?

mu_direction = 0                #mu for von mises distribution selecting direction of movements

kappa_direction = 0             #sigma for von mises distribution

mu_distance = 4                 #mean movement-distance (lognormal distribution)

sigma_distance = 0.1            #sd of movement distance

sex = True                      #is this a sexual species?

repro_age = (8, 5)           #age at sexual maturity (int or float for non-sexual species, tuple or list of two ints/floats for sexual species; set to 'None' to not make this an age-structured species


mating_radius = 50              #radius of mate-searching area

mu_dispersal = 5                #mean dispersal distance (lognormal distribution)

sigma_dispersal = 0.4          #sd of dispersal distance


size = [1] * T                  # float/int, or list/tuple of length T containing floats/ints, expressingthe target population size over model time as a ratio of the starting size (N)


#########################
      ###################
      ###################
#########################


#set seed, if requested
if set_seed:
    random.seed(seed_num)
    r.seed(seed_num)





#if run in debugging mode, set debug = true
debug = False
if '-d' in sys.argv or '--debug' in sys.argv:
    debug = True


#if running plots requested, set running_plot = true

running_plot = False
if '--plots' in sys.argv:
    running_plot = True




#read parameter file
param_file = "???" #replace with command-line arg
params = readcsv(param_file)

#parse params into appropriate structure






#------------------------------------------------#                         


#---------------------#
# ESTABLISH LANDSCAPE #
#---------------------#



#if landscape/s to be pulled from file, do so




#if landscape/s to be generated, then generate it/them
if '-l' in sys.argv or '--landscape' in sys.argv: #OR SOMETHING IN PARAM FILE

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



pop = population.create_population(N = N, genomic_arch = genomic_arch, dims = dims, size = size, T = T) 





#pickle initial pop, for later analysis
pop.pickle('initial_pop.p')


#------------------------------------------------#                         
                                                                                                                                                               
#---------#
# BURN IN #
#---------#


for t in range(1, burn_T):
    if move:

        [movement.move(ind, land, mu_direction = mu_direction, kappa_direction = kappa_direction, mu_distance = mu_distance, sigma_distance = sigma_distance) for ind in pop.individs.values()];
        

    
    #find mates
    mating_pairs = pop.find_mates(mating_radius)


    #mate

    #cross-over and produce gametes

    for pair in mating_pairs:
        #create zygote from the mating event
        zygote = mating.mate(pop, pair, genomic_arch)


        #generate new offspring


        #disperse offspring
    
    #mutate
    
    #select

    #add to ages
    pop.birthday()



#pickle pop after burn-in, for later analysis
pop.pickle('burnt_pop.p')


#------------------------------------------------#                         
                                                                                                                                                               
#-----------#
# RUN MODEL #
#-----------#

#once the starting, burned-in population is created, run the model for each of N different parameter sets in the input file


for N in param_sets:

    for t in range(1, T):
        if move:
            pass
            # move

            #find mates
            #mate
            #cross-over
            #disperse offspring
            #mutate
            #select




#------------------------------------------------#
                                                                                                                                                               
#-----------------------------#
# POST-PROCESS RESULTING DATA #
#-----------------------------#

#post-process resulting data, to output in desired format (as specified in command-line/param-file args

