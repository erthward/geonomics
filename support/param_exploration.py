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
import matplotlib.pyplot as plt
from copy import deepcopy

from ops import movement



#------------------------------------
# CLASSES ---------------------------
#------------------------------------


#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------



def plot_movement(pop, land, scape_num, num_timesteps, params = None, mu_direction = None, kappa_direction =
        None, mu_distance = None, sigma_distance = None, move_surf = None, subset_pop = None, color = 'black', color_by_individ = False):

    '''Useful for visual exploration of the movement parameters. Either provide a params dictionary with the
    necessary information in it, or provide the necessary parameters piecemeal.'''

    if params == None and mu_distance == None:
        print('\n\n\tERROR: Must provide either a params dictionary or else stipulate parameter values.\n\tPlease try again.\n\n')
        return None

    #if a params dictionary was fed in, copy it as toy_params
    if params is not None:
        toy_params = deepcopy(params)


        toy_m_params = toy_params.comm.pops[0].movement
        #then override existing values in params dictionary, if a new value has been entered for trial
        if mu_direction != None:
            toy_m_params['mu_direction'] = mu_direction
        if kappa_direction != None:
            toy_m_params['kappa_direction'] = kappa_direction
        if mu_distance  != None:
            toy_m_params['mu_distance'] = mu_distance 
        if sigma_distance != None:
            toy_m_params['sigma_distance'] = sigma_distance
        if move_surf != None:
            toy_m_params['move_surf'] = move_surf

        toy_params.comm.pops[0].movement = toy_m_params


    #if no params dictionary provided, create toy_params from the individual parameter values that were fed in
    elif params is None:
        toy_m_params = {}
        toy_m_params['mu_direction'] = mu_direction
        toy_m_params['kappa_direction'] = kappa_direction
        toy_m_params['mu_distance'] = mu_distance 
        toy_m_params['sigma_distance'] = sigma_distance
        if move_surf != None:
            toy_m_params['move_surf'] = move_surf
        elif land.move_surf == None:
            toy_m_params['move_surf'] = False
        else:
            toy_m_params['move_surf'] = True

        toy_params = {'pop':{'movement':toy_m_params}}



    toy_pop = deepcopy(pop)


    if subset_pop:
        cull_individs = np.random.choice([*toy_pop], len(toy_pop) - subset_pop, replace = False)
        [toy_pop.pop(ind) for ind in cull_individs]


    toy_pop.plot(land, scape_num, color = 'white', size = 10)

    linewidths = np.linspace(2,5, num = num_timesteps)

    #Create a colors list, to plot different-color paths for different individuals
    #colors = ['black', 'red', 'orange', 'yellow', 'green', 'blue', 'white']
    colors = [plt.cm.Accent(i) for i in np.linspace(0, 0.9, 9)]


    #NOTE: offset all values by -0.5 for visual reconciliation, to account for apparent offset of axes on top of raster?
    new_x = [ind.x-0.5 for ind in toy_pop.values()]
    new_y = [ind.y-0.5 for ind in toy_pop.values()]


    for t in range(num_timesteps):
        old_x = [x for x in new_x]
        old_y = [y for y in new_y]

        #NOTE: offset all values by -0.5 for visual reconciliation, to account for apparent offset of axes on top of raster?
        toy_pop.move(land, toy_params)
        new_x = [ind.x-0.5 for ind in toy_pop.values()]
        new_y = [ind.y-0.5 for ind in toy_pop.values()]

        if color_by_individ == True:
            [mpl.pyplot.plot((old_x[i], new_x[i]), (old_y[i], new_y[i]), '-', scalex = False, scaley = False, 
            linewidth = linewidths[t], color = colors[i%len(colors)], alpha = 0.5) for i in range(len(old_x))];

        else:
            mpl.pyplot.plot((old_x, new_x), (old_y, new_y), '-', scalex = False, scaley = False, 
            linewidth = linewidths[t], color = color, alpha = 0.5)

