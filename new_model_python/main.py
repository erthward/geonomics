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

import numpy as np
from numpy import random as r
import random


import genome
import population




#------------------------------------------------#                         


#-----------------#
# GRAB PARAMETERS #
#-----------------#


#grab (and process some) command-line arguments, setting majority of "secondary" params to their default
#values ("typical", "characteristic", or uniformly-distributed values) unless otherwise stipulated in the
#command-line args



#set seed, if requested
if set_seed:
    random.seed(1)
    r.seed(1)





#########################
      ###################
# TEMP PARAMS FOR TESTING
      ###################
#########################

T = 100                         #total model runtime

L = 1e5                         #total number of loci

n = 23                          #number of chromosomes

N = 700                         #total pop size

dims = (100,100)      #dimensions of landscape  


#########################
      ###################
      ###################
#########################





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



pop = population.create_population(N, genomic_arch, dims, 0.2, 0.2) 




#------------------------------------------------#                         
                                                                                                                                                               
#---------#
# BURN IN #
#---------#

for t in range(1, burn_T)
  if move
    # move
  end
  if n%gen_t == 0
    #find mates
    #mate
    #cross-over
    #disperse offspring
    #mutate
    #select
  end
end



#------------------------------------------------#                         
                                                                                                                                                               
#-----------#
# RUN MODEL #
#-----------#

#once the starting, burned-in population is created, run the model for each of N different parameter sets in the input file


for N in param_sets

  for t in range(1, T)
    if move
      # move
    end
    if n%gen_t == 0
      #find mates
      #mate
      #cross-over
      #disperse offspring
      #mutate
      #select
    end
  end

end




#------------------------------------------------#
                                                                                                                                                               
#-----------------------------#
# POST-PROCESS RESULTING DATA #
#-----------------------------#

#post-process resulting data, to output in desired format (as specified in command-line/param-file args

