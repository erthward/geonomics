#!/usr/bin/python
#param_exploration.py


'''
##########################################

Module name:                param_exploration.py

Module contains:
                            - functions to assist in exploring key parameters


Author:                     Drew Ellison Hart
Email:                      drew.hart@berkeley.edu
Github:                     URL
Start date:                 12-28-15
Documentation:              URL


##########################################
'''


import numpy as np
from numpy.random import vonmises, lognormal
import matplotlib as mpl

import movement



#------------------------------------
# CLASSES ---------------------------
#------------------------------------


#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------



def watch_movement(pop, land, scape_num, mu_direction, kappa_direction, mu_distance, sigma_distance, num_timesteps, resist_surf = None):

    pop.show(land, scape_num)

    new_x = [ind.x for ind in pop.individs.values()]
    new_y = [ind.y for ind in pop.individs.values()]

    for i in range(num_timesteps):
        old_x = [x for x in new_x]
        old_y = [y for y in new_y]
        new_coords = [movement.move(ind, land, mu_direction, kappa_direction, mu_distance, sigma_distance, resist_surf) for ind in pop.individs.values()]
        new_x = [coord[0] for coord in new_coords]
        new_y = [coord[1] for coord in new_coords]

        mpl.pyplot.plot((old_x, new_x), (old_y, new_y), 'k-', scalex = False, scaley = False, color = 'black')





   

