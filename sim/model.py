#!/usr/bin/python
# model.py

'''
##########################################

Module name:          sim.model


Module contains:
                      - classes and functions to facilitate building a Geonomics model and 
                      running multiple iterations of it, including across multiple cores


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           07-06-18
Documentation:        URL


##########################################
'''

#geonomics imports
from structs import landscape, community, genome
from sim import burnin

#other imports
import numpy.random as r
import random
import copy


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################



#TODO:
    # finish and debug functions I've written
    # add verbose statements
    # finish the data and stats modules and tie them in correctly


#the Model class, which will contain the main model objects, build the
#burn-in and main functions, allow customization of the objects, manage and run
#iterations, and produce and save required stats and data
class Model:
    def __init__(self, params, verbose=False):

        #get the model params
        m_params = params.model

        #set the seed, if required
        self.seed = None
        if 'seed' in [*m_params]:
            if m_params.seed.set:
                self.seed = m_params.seed.num
            self.set_seeds()

        #get minimum burn-in runtime (in timesteps)
        self.burn_T = m_params.time.burn_T
        #and set the burn-timestep counter (burn_t) to 0
        self.burn_t = 0
        #get main runtime (in timesteps)
        self.T = params.model.time.T
        #and set the timestep counter (t) to 0
        self.t = 0

        #make the land and community objects
        self.land = self.make_land(params)
        self.comm = self.make_community(self.land, params)

        #get the number of model iterations to run
        self.n_its = m_params.its.n_its
        #set the its list
        self.its = [*range(self.n_its)][::-1]
        #start the it counter
        self.it = None
        #and set it
        self.set_it()

        #and set the orig_land and orig_comm objects, if rand_land and rand_comm are False
        self.rand_land = m_params.its.rand_land
        self.rand_comm = m_params.its.rand_comm
        self.orig_land = None
        self.orig_comm = None
        if not self.rand_land:
            self.orig_land = deepcopy(self.land)
        if not self.rand_comm:
            self.orig_comm = deepcopy(self.comm)

        #grab the data params and stats params
        self.data_params = copy.deepcopy(m_params.data)
        self.stats_params = copy.deepcopy(m_params.stats)

    #method to set the curr_it counter
    def set_it(self):
        self.it = self.n_its.pop()

    #method to set seed (will be run when Model object is first created, if
    #called for in params)
    def set_seeds(self):
        random.seed(self.seed)
        r.seed(self.seed)

    #method to wrap around landscape.make_land
    def make_land(self, params):
        landscape.make_land(params)

    #method to wrap around community.make_community
    def make_community(self, params):
        community.make_community(self.land, params, burn = True)

    #method to create a sim function queue (main sim fn, unless burn is True)
    def make_sim_fn_queue(self, burn=False, reassign_genomes=True):
    #NOTE: the queue will be built as a list of functions, so 
    #that the user can add, remove, change the order, or repeat functions
    #within the list as desired
    queue = []
    #reassign the genomes if it's a main function and reassign_genomes is True (default)
    if not burn:
        if reassign_genomes:
            #createa a lambda fn with genome.set_genomes, so that the number of
            #mutations will be properly estimated at the time that the function
            #is called, after burn-in births have happened
            for pop in self.comm.values():
                reassign_fn = lambda self: genome.set_genomes(pop, self.burn_T, self.T)
                queue.append(reassign_fn) 
    #append the set_age_stage methods
    for pop in self.comm.values():
        queue.append(pop.set_age_stage)
    #append the set_Nt methods
    for pop in self.comm.values():
        queue.append(pop.set_Nt)
    #append the do_movement_methods, if pop.move
    for pop in self.comm.values():
        if pop.move:
            queue.append(lambda: pop.do_movement(land))
    #append the do_pop_dynamics methods
    #FIXME: Consider whether the order of these needs to be specified, or
    #randomized, should people want to eventually simulate multiple, interacting populations
    for pop in self.comm.values():
        queue.append(lambda: pop.do_pop_dynamics(land, with_selection = (pop.selection and not burn)))

    #add the make_change methods, if not burn-in and if pop and/or land have Changer objects
    if not burn:
        #add land.changer method
        if land.changer is not None:
            queue.append(lambda: land.make_change(self.t))
        for pop in self.comm.values():
            if pop.changer is not None:
                queue.append(pop.make_change)

    #TODO add the data.save and stats.save methods, if not burn-in and if needed
            
    #TODO add the burn-in function if need be
    if burn:
        burn_fn = None
        queue.append(burn_fn)
    
  
    #TODO a method to set the model's next iteration
    def set_next_iteration(self, core=None):
        #deepcopy the original land if not None, else randomly generate new land
        
        #deepcopy the original comm if not None, else randomly generate new land

        #update the iteration numbers

        #reset the self.burn_t and self.t attributes to 0

        #create burn-in fns, if this iteration's community object wasn't
        #deepcopied, or WAS but has burnt = False

        #create Data and Stats objects, if required, using self.data_params and
        #self.stats_params
            #and call their methods to create new data and stats subdirectories for
            #this iteration
    
    
    #TODO a method for running the model's next iteration
    def run_next_iteration(self):
        #set the next iteration to be run
        self.set_next_iteration()

        #loop over the burn-in timesteps running the burnin queue (if copied pop isn't already burned in)

        #set the pop.t attribute to 0, and pop.T attribute to self.T

        #loop over the timesteps, running the main function repeatedly
            #drop the first function after timestep 0, if it's the reassign genomes fn

            #end the iteration if any population's pop.t attribute is -1 (which indiates extinction)

    #TODO a method to run the overall model
    def run(self, verbose=False):
        pass


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

