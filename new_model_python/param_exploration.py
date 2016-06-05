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



def watch_movement(pop, land, scape_num, num_timesteps, params = None, mu_direction = None, kappa_direction = None, mu_distance = None, sigma_distance = None, resist_surf = None, subset_pop = None):

    if params:
        mu_direction = params['mu_direction']
        kappa_direction = params['kappa_direction']
        mu_distance = params['mu_distance']
        sigma_distance = params['sigma_distance']

    import copy

    toy_pop = copy.deepcopy(pop)

    if subset_pop:
        cull_individs = np.random.choice(toy_pop.individs.keys(), len(toy_pop.individs) - subset_pop, replace = False)
        [toy_pop.individs.pop(ind) for ind in cull_individs]


    toy_pop.show(land, scape_num)

    new_x = [ind.x for ind in toy_pop.individs.values()]
    new_y = [ind.y for ind in toy_pop.individs.values()]

    for i in range(num_timesteps):
        old_x = [x for x in new_x]
        old_y = [y for y in new_y]
        [ind.move(land, mu_direction, kappa_direction, mu_distance, sigma_distance, resist_surf) for ind in toy_pop.individs.values()];
        new_x = [ind.x for ind in toy_pop.individs.values()]
        new_y = [ind.y for ind in toy_pop.individs.values()]

        mpl.pyplot.plot((old_x, new_x), (old_y, new_y), 'k-', scalex = False, scaley = False, color = 'black')





   

