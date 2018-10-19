#!/usr/bin/python 
# community.py

'''
##########################################

Module name:              community


Module contains:
                          - definition of the Community class
                          - a make_commmunity() function
                          - nothing else for now, but it leaves open the possibility of adding functionality for mutliple populations (e.g. species interactions, speciation models, etc.)


Author:                    Drew Ellison Hart
Email:                     drew.hart@berkeley.edu
Github:                    URL
Start date:                07-25-18
Documentation:             URL


##########################################
'''

#geonomics imports
from structs import population
from sim import burnin

#other imports
import numpy as np


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

#NOTE: this class won't be doing too much right away, but it lays the foundation 
#for enabling interactions between populations (i.e. species) further down the road
class Community(dict):
    def __init__(self, land, pops):
        self.land = land
        self.update(pops)
        self.n_pops = len(pops)
        self.t = -1 #counter for timesteps (starts at -1, to indicate that the
                    #community is unrun, and so that the first timestep will be set to 0 at
                    #beginning of the timestep)
        #set the burned attribute (defaults to False, but will be set to True after burn-in 
        #has successfully completed, and will be used by the Model object to determine whether 
        #or not burn-in needs to happen each iteration)
        self.burned = False

    #define the __str__ and __repr__ special methods
    def __str__(self):
        #get the longest pop name and longest string-repr pop size, to be used
        #to horizontally rectify all the lines for each of the pops
        max_len_pop_name = max([len(pop.name) for pop in self.values()])
        max_len_pop_size = max([len(str(len(pop))) for pop in self.values()])
        #get a string representation of the class
        type_str = str(type(self))
        #get a string representation of the populations
        pops_str = '%i Population%s:\n' % (len(self), 's' * (len(self) > 1))
        pops_str = pops_str + '\n'.join(['\tpop %i: ' % k +
            "%s'%s' (%s%s inds.): " % (' ' * (max_len_pop_name - len(v.name)),
            v.name, ' ' * (max_len_pop_size - len(str(len(v)))),
            str(len(v))) +
            str(v).split('\n')[0] for k,v in self.items()])
        ##get a string representation of the first two and last two parameters
        #params_str = "\nParameters:\n\t" + ',\n\t'.join(sorted([str(k) +
        #    ': ' + str(v) for k,v in vars(self).items()][:2])) + ','
        #params_str = params_str + '\n\t...\n\t'
        #params_str = params_str + ',\n\t'.join(sorted([str(k) +
        #    ': ' + str(v) for k,v in vars(self).items()][-2:]))

        #return '\n'.join([type_str, pops_str, params_str])
        return '\n'.join([type_str, pops_str])

    def __repr__(self):
        repr_str = self.__str__()
        return repr_str



    #method to increment the self.t attribute (the timestep counter)
    def _set_t(self):
        self.t += 1

    #method to reset the self.t attribute (the timestep counter)
    def _reset_t(self):
        self.t = -1

    #method to check if all populations have burned in
    def _check_burned(self, burn_T):
        #check minimum burn-in time has passed
        burnin_status = np.all([len(pop.Nt) >= burn_T for pop in self.values()])
        #if so, then check all burn-in tests for all pops
        if burnin_status:
            adf_tests = np.all([burnin._test_adf_threshold(pop, burn_T) for
                                pop in self.values()])
            t_tests = np.all([burnin._test_t_threshold(pop, burn_T) for
                                pop in self.values()])
            burnin_status = adf_tests and t_tests
        #set the community burn-in tracker
        self.burned = burnin_status
        #and set the pops' burn-in trackers
        [setattr(pop, 'burned', burnin_status) for pop in self.values()]


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#function for making a community using a Landscape object and params
def _make_community(land, params, burn=False):
    pops = {n: population._make_population(land = land, name = name, 
        pop_params = params.comm.pops[name], burn = burn) for n, name in
            enumerate(params.comm.pops.keys())}
    return Community(land, pops) 

