#!/usr/bin/python
#movement.py


'''
##########################################

Module name:                movement

Module contains:
                            - function for simulating the movement of an individual, according to input parameters
                            - associated functions


Author:                     Drew Ellison Hart
Email:                      drew.hart@berkeley.edu
Github:                     URL
Start date:                 12-28-15
Documentation:              URL


##########################################
'''

import numpy as np
from numpy.random import vonmises, lognormal


#----------------------------------
# TYPES ---------------------------
#----------------------------------






#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------



def move(individual, land, mu_direction, kappa_direction, mu_distance, sigma_distance, resist_surf = None):
    if resist_surf:
        pass #NOTE: this should be a function that allows individuals' movements to be determined by some calculation on a resistance surface raster
    #yet to figure this out; was thinking perhaps of surrounding the surface with a ring of cells with value 0, then calculating a Von Mises mixture distribution (a mixture of all 8 queen's neighbourhood cells) for each cell, then using that in each cell to choose the direction of an individual's movement...
    else:
        #Implement a random walk
        direction = vonmises(mu_direction, kappa_direction)
        #NOTE: Instead of lognormal, could use something with long right tail for Levy-flight type movement
        distance = lognormal(mu_distance, sigma_distance) 
        new_x = min(max(individual.x + np.cos(direction)*distance, 0), land.dims[0]-0.001)
        new_y = min(max(individual.y + np.sin(direction)*distance, 0), land.dims[1]-0.001)
        
        return new_x, new_y




