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

        #get the full input params dict
        self.params = deepcopy(params)
        #get the model params (to make following code shorter)
        m_params = self.params.model

        #grab the data params and stats params
        #self.data_params = copy.deepcopy(m_params.data)
        #self.stats_params = copy.deepcopy(m_params.stats)

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

        #get the number of model iterations to run
        self.n_its = m_params.its.n_its
        #set the its list
        self.its = [*range(self.n_its)][::-1]
        #start the it counter
        self.it = None
        #and set it
        self.set_it()

        #make the land and community objects
        self.land = self.make_land(params)
        self.comm = self.make_community(params)

        #create a self.reassign_genomes attribute, which defaults to False,
        #unless any population has a genomic architecture, indicating that its
        #genomes should be reassigned after burn-in; in that case it will be 
        #reset to False as soon as the genomes are reassigned
        self.reassign_genomes = np.any([pop.gen_arch is not None for pop in self.comm.pops])

        #and set the orig_land and orig_comm objects, if rand_land and rand_comm are False
        self.rand_land = m_params.its.rand_land
        self.rand_comm = m_params.its.rand_comm
        self.orig_land = None
        self.orig_comm = None
        if not self.rand_land:
            self.orig_land = deepcopy(self.land)
        if not self.rand_comm:
            self.orig_comm = deepcopy(self.comm)

        #and set the self.rand_burnin attribute; if this is False and
        #self.rand_comm is False, self.orig_comm will be replaced with the
        #first iteration's community as soon as it has exited burn-in and had
        #its genomes reassigned
        self.rand_burnin = m_params.its.rand_burnin

        #create the burn-in and main queues
        self.burnin_fn_queue = self.make_fn_queue(burn = True)
        self.main_fn_queue = self.make_fn_queue(burn = False)


    #method to set the iteration-counter, self.it
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

    #method to create the simulation functionality, as a function queue 
    #NOTE: (creates the list of functions that will be run by
    #self.run_burnin_timestep or self.run_main_timestep, depending on burn argument)
    #NOTE: because the queue is a list of functions, the user can add, remove, 
    #change the order, or repeat functions within the list as desired
    def make_fn_queue(self, burn=False):
        queue = []
        #append the set_age_stage methods to the queue
        for pop in self.comm.pops:
            queue.append(pop.set_age_stage)
        #append the set_Nt methods
        for pop in self.comm.pops:
            queue.append(pop.set_Nt)
        #append the do_movement_methods, if pop.move
        for pop in self.comm.pops:
            if pop.move:
                queue.append(pop.do_movement)
        #append the do_pop_dynamics methods
        #FIXME: Consider whether the order of these needs to be specified, or
        #randomized, should people want to eventually simulate multiple, interacting populations
        for pop in self.comm.pops:
            queue.append(pop.do_pop_dynamics)

        #add the make_change methods, if not burn-in and if pop and/or land have Changer objects
        if not burn:
            #add land.changer method
            if self.land.changer is not None:
                queue.append(lambda: self.land.make_change(self.t))
            for pop in self.comm.pops:
                if pop.changer is not None:
                    queue.append(pop.make_change)


        #TODO add the data.save and stats.save methods, if not burn-in and if needed


        #add the burn-in function if need be
        if burn:
            for pop in self.comm.pops:
                queue.append(lambda: pop.check_burned(self.burn_T))

        return(queue)


    #a method to run the main function queue (and reassign genomes, if necessary)
    def run_main_timestep(self):
        #reassign the genomes, if self.reassign_genomes indicates
        if self.reassign_genomes:
            for pop in self.comm.pops:
                if pop.gen_arch is not None:
                    genome.set_genomes(pop, self.burn_T, self.T)
            #and then set the model.reassign_genomes attribute to False, so
            #that they won'r get reassigned again during this iteration
            model.reassign_genomes = False
        #call all of the functions in the main function queue
        [fn() for fn in self.main_fn_queue]


    #define the burn-in function
    def run_burnin_timestep(self):
        #call all of the functions in the burn-in function queue
        [fn() for fn in self.burnin_fn_queue]


    #TODO a method to set the model's next iteration
    def set_next_iteration(self, core=None):
        #deepcopy the original land if not None, else randomly generate new land

        #deepcopy the original comm if not None, else randomly generate new land

        #update the iteration numbers

        #reset the self.burn_t and self.t attributes to 0

        #create burn-in fns, if this iteration's community object wasn't
        #deepcopied, or WAS but has burned = False

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

        #loop over the timesteps, running the run_main function repeatedly
            #end the iteration if any population's pop.t attribute is -1 (which indiates extinction)


    #TODO a method to run the overall model
    def run_model(self, verbose=False):
        pass
        #TODO: Add some handler for farming instances out to different nodes, if available?

        #for n in range(self.n_its):

            #if verbose:
                #print iteration number

            #set the next iteration

            #check if self.orig_comm and self.orig_land, then either randomly
            #create new ones or just deepcopy those for the next iteration, as appropriate

            #while not np.all([pop.burned for pop in self.comm.pops]): 
                #self.burn_t += 1
                #if verbose:
                    #print that it's burning in, and print some output at each timestep
                #burn in the next iteration


            #if self.orig_comm is not None and not self.rand_burn:
                #switch out self.orig_comm for the now burned-in comm

            #for t in range(self.T):
                #if verbose:
                    #print that it's now running the main model, and print some
                    #output at each timestep
                #self.t = t 
                #and set self.comm.t and self.comm.pop.t attributes as well?
                #run the main_queue
                #check extinct, break and report if so


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################


