#!/usr/bin/python
# movement.py


'''
##########################################

Module name:                movement

Module contains:
                            - function for simulating the movement of an
                              individual, according to input parameters
                            - associated functions


Author:                     Drew Ellison Hart
Email:                      drew.hart@berkeley.edu
Github:                     URL
Start date:                 12-28-15
Documentation:              URL


##########################################
'''


#TODO:
#   - vectorize dispersal (i.e. create all offspring and parent midpoints,
#     then draw new locations for all offspring simultaneously)
#   - create simpler (private?) methods for making directional and distance
#     draws, then call those within move, disperse, etc



import numpy as np
import numpy.random as r
import random
from numpy import sin, cos, pi
from numpy.random import vonmises as r_vonmises
from numpy.random import lognormal, gamma, wald
import matplotlib.pyplot as plt
from scipy.stats import vonmises as s_vonmises
from collections import Counter as C

s_vonmises.a = -np.inf
s_vonmises.b = np.inf


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

def _move(spp):
    #get individuals' coordinates (soon to be their old coords, so 
    #'old_x' and 'old_y')
    old_x, old_y = [a.flatten() for a in np.split(spp._get_coords(),
                                                            2, axis = 1)]
    #and get their cells (by rounding down to the int)
    old_x_cells, old_y_cells = [a.flatten() for a in np.split(spp._get_cells(),
                                                              2, axis = 1)]
    # choose direction using movement surface, if applicable
    if spp._move_surf:
        #and use those choices to draw movement directions
        direction = spp._move_surf._draw_directions(old_x_cells, old_y_cells)
        # NOTE: Pretty sure that I don't need to constrain values output
        #for the Gaussian KDE that is approximating the von Mises mixture 
        #distribution to 0<=val<=2*pi, because e.g. cos(2*pi + 1) = cos(1),
        #etc...
        # NOTE: indexed out of move_surf as y then x because becuase the 
        #list of lists (like a numpy array structure) is indexed i then j,
        #i.e. vertical, then horizontal
    # else, choose direction using a random walk with a uniform vonmises
    elif not spp._move_surf:
        direction = r_vonmises(spp.direction_distr_mu,
                               spp.direction_distr_kappa, size = len(old_x))

    # choose distance
    # NOTE: Instead of lognormal, could use something with long right tail
    #for Levy-flight type movement, same as below
    distance = wald(spp.distance_distr_mu, spp.distance_distr_sigma,
                                                        size = len(old_x))

    #create the new locations by adding x- and y-dim line segments to their
    #current positions, using trig then clip the values to be within the 
    #landscape dimensions
    #NOTE: subtract a small value to avoid having the dimension itself set
    #as a coordinate, when the coordinates are converted to np.float32 
    new_x = old_x + cos(direction)*distance
    new_x = np.clip(new_x, a_min = 0, a_max = spp._land_dim[1]-0.001)
    new_y = old_y + sin(direction)*distance
    new_y = np.clip(new_y, a_min = 0, a_max = spp._land_dim[0]-0.001)
    #then feed the new locations into each individual's set_pos method
    [ind._set_pos(x, y) for ind, x, y in zip(spp.values(), new_x, new_y)];


def _disperse(spp, land, parent_midpoint_x, parent_midpoint_y,
        dispersal_distr_mu, dispersal_distr_sigma, mu_dir = 0, kappa_dir = 0):
    within_landscape = False
    while within_landscape == False:
        # choose direction using movement surface, if applicable
        if spp._disp_surf:
            #and use those choices to draw movement directions
            direction = spp._disp_surf._draw_directions(
                [int(parent_midpoint_x)], [int(parent_midpoint_y)])[0]
        # else, choose direction using a random walk with a uniform vonmises
        elif not spp._disp_surf:
            direction = r_vonmises(mu_dir, kappa_dir)
        distance = wald(dispersal_distr_mu, dispersal_distr_sigma)
        offspring_x = parent_midpoint_x + np.cos(direction)*distance
        offspring_y = parent_midpoint_y + np.sin(direction)*distance
        offspring_x = np.clip(offspring_x, a_min =0, a_max = land.dim[1]-0.001)
        offspring_y = np.clip(offspring_y, a_min =0, a_max = land.dim[0]-0.001)
        within_landscape = (offspring_x > 0
                            and offspring_x < land.dim[1]) and (offspring_y > 0
                                                and offspring_y < land.dim[0])
    return (offspring_x, offspring_y)

