    @@@@@  @@@@@  @@@@@       @    @  @@@@@  @     @  @@@@@  @@@@@  @@@@@
   @      @      @   @       @@   @  @   @  @@   @@    @    @      @    
  @  @@  @@@@@  @   @  @@@@ @ @  @  @   @  @ @ @ @    @    @      @@@@@
 @   @  @      @   @       @  @ @  @   @  @  @  @    @    @          @
@@@@@  @@@@@  @@@@@       @    @  @@@@@  @     @  @@@@@  @@@@@  @@@@@


#=
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
  

=#




#if -h/--help flag in command-line args, return docstring and exit


#------------------------------------------------#



#---------------#
# LOAD PACKAGES #
#---------------#

using Distributions
using genome





#------------------------------------------------#                         


#-----------------#
# GRAB PARAMETERS #
#-----------------#


#grab (and process some) command-line arguments, setting majority of "secondary" params to their default
#values ("typical", "characteristic", or uniformly-distributed values) unless otherwise stipulated in the
#command-line args



#set seed, if requested
if set_seed 
  srand(1)
end





#########################
      ###################
# TEMP PARAMS FOR TESTING
      ###################
#########################

L = 1e4  #total number of loci

n = 5    #number of chromosomes

N = 700  #total pop size

#########################
      ###################
      ###################
#########################





#if run in debugging mode, set debug = true


#if running plots requested, set running_plot = true





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
using landscape_creation

landscapes = Dict()

for i in range(1, num_landscapes)
  landscapes[i] = gen_land(params) 
end



  
#------------------------------------------------#                         
                                                                                                                                                               
#-----------------------------#
# CREATE GENOMIC ARCHITECTURE #
#-----------------------------#

using genome

genomic_arch, dist_s, dist_D = genome.build_genomic_arch(n,l_c)



                                                                                                                                                               
#------------------------------------------------#



#----------------------#
# ESTABLISH POPULATION #
#----------------------#

using population


pop = population.sim_population(N) #




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

