#!/usr/bin/python
#param_help.py


'''
##########################################

Module name:                param_help.py

Module contains:
                            - Functions to assist in exploring key parameters


Author:                     Drew Ellison Hart
Email:                      drew.hart@berkeley.edu
Github:                     URL
Start date:                 12-28-15
Documentation:              URL


##########################################
'''

#geonomics imports
from ops import movement

#other imports
import numpy as np
from numpy.random import vonmises, lognormal
import matplotlib as mpl
import matplotlib.pyplot as plt
from copy import deepcopy

#------------------------------------
# CLASSES ---------------------------
#------------------------------------


#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------



def plot_movement(spp, land, num_timesteps, lyr_num=None, params=None,
        direction_distr_mu=None, direction_distr_kappa=None,
        distance_distr_mu=None, distance_distr_sigma=None, move_surf=None,
        subset_spp=None, color='black', color_by_individ=False, size=10):

    '''
    Useful for visual exploration of the movement parameters.
    Either provide a params dictionary with the necessary
    information in it, or provide the necessary parameters piecemeal.
    '''

    #assert that either a params dict was provided or parameters were
    #stipulated
    assert (params is not None
        or (distance_distr_mu is not None
            and distance_distr_sigma is not None
            and direction_distr_mu is not None
            and direction_distr_kappa is not None)), ("Either a parameters "
        "dictionary must be provided, or parameters must be individually "
        "defined as arguments.")

    #if a params dictionary was fed in, copy it as toy_params
    if params is not None:
        toy_params = deepcopy(params)
        #copy the movement params from toy_params
        toy_move_params = toy_params.comm.species[spp.name].movement
        #then override existing values in the movement params dictionary,
        #if they have new values that have been entered for this trial
        if direction_distr_mu != None:
            toy_move_params['direction_distr_mu'] = direction_distr_mu
        if direction_distr_kappa != None:
            toy_move_params['direction_distr_kappa'] = direction_distr_kappa
        if distance_distr_mu  != None:
            toy_move_params['distance_distr_mu'] = distance_distr_mu 
        if distance_distr_sigma != None:
            toy_move_params['distance_distr_sigma'] = distance_distr_sigma
        if move_surf != None:
            toy_move_params['move_surf'] = move_surf
        toy_params.comm.species[spp.name].movement = toy_move_params

    #if no params dictionary provided, create toy_params from
    #the individual parameter values that were fed in
    elif params is None:
        toy_move_params = {}
        toy_move_params['direction_distr_mu'] = direction_distr_mu
        toy_move_params['direction_distr_kappa'] = direction_distr_kappa
        toy_move_params['distance_distr_mu'] = distance_distr_mu 
        toy_move_params['distance_distr_sigma'] = distance_distr_sigma
        if move_surf != None:
            toy_move_params['move_surf'] = move_surf
        elif spp.move_surf == None:
            toy_move_params['move_surf'] = False
        else:
            toy_move_params['move_surf'] = True
            toy_params = {'comm':{'species':{'movement':toy_move_params}}}

    #now copy the species to a toy object
    toy_spp = deepcopy(spp)
    #replace its movement-parameter values
    for param in ['direction_distr_mu', 'direction_distr_kappa',
        'distance_distr_mu', 'distance_distr_sigma']:
        setattr(toy_spp, param,
            toy_params['comm']['species'][spp.name]['movement'][param])

    #replace its movement surface, if necessary
    if move_surf is not None:
        toy_spp._move_surf = move_surf

    #subset individuals, if necessary
    if subset_spp:
        cull_individs = np.random.choice([*toy_spp],
            len(toy_spp) - subset_spp, replace = False)
        [toy_spp.pop(ind) for ind in cull_individs]
        #reset the coords before plotting
        toy_spp._set_coords_and_cells()

    #set the plotting linewidths to increase over runtime
    linewidths = np.linspace(2,5, num = num_timesteps)

    #Create a colors list, to plot different-colored paths for individuals
    #colors = ['black', 'red', 'orange', 'yellow', 'green', 'blue', 'white']
    colors = [plt.cm.Accent(i) for i in np.linspace(0, 0.9, 9)]

    #plot the species at its starting locations
    toy_spp._plot(lyr_num=lyr_num, land = land, color = 'white', size = 10)

    #set the new_x and new_y objects to the current locations, before for loop
    #that will iteratively update them to the newer locations after movement
    new_x = [ind.x for ind in toy_spp.values()]
    new_y = [ind.y for ind in toy_spp.values()]

    #loop for the number of timesteps, iteratively moving individuals and
    #plotting the vectors between their previous and new positions
    for t in range(num_timesteps):
        old_x = [x for x in new_x]
        old_y = [y for y in new_y]

        toy_spp._do_movement(land)

        new_x = [ind.x for ind in toy_spp.values()]
        new_y = [ind.y for ind in toy_spp.values()]

        #plot with separate colors for each individual, if necessary
        if color_by_individ == True:
            [plt.plot((old_x[i], new_x[i]), (old_y[i], new_y[i]), '-',
                scalex = False, scaley = False, linewidth = linewidths[t],
                color = colors[i%len(colors)],
                alpha = 0.5) for i in range(len(old_x))];

        #or else plot without individuals colored differently
        else:
            mpl.pyplot.plot((old_x, new_x), (old_y, new_y), '-',
                scalex = False, scaley = False, linewidth = linewidths[t],
                color = color, alpha = 0.5)

