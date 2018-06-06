#!/usr/bin/python
# movement.py


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
import numpy.random as r
import random
from numpy import sin, cos, pi
from numpy.random import vonmises as r_vonmises
from numpy.random import lognormal, choice, gamma, wald
import matplotlib.pyplot as plt
from scipy.stats import vonmises as s_vonmises
import copy
from operator import itemgetter
from numba import jit
from collections import Counter as C

s_vonmises.a = -np.inf
s_vonmises.b = np.inf


# ----------------------------------
# TYPES ---------------------------
# ----------------------------------


# --------------------------------------
# FUNCTIONS ---------------------------
# --------------------------------------

def move(pop, land, params):
    
    #get individuals' ids
    #ids = pop.individs.keys()
    #plug them into an itemgetter
    #ind_ig = itemgetter(*ids)
    #then use the itemgetter to get individuals' coordinates
    old_x, old_y = [a.flatten() for a in np.split(pop.get_coords(), 2, axis = 1)]
    #and their cells (by rounding down to the int)
    old_x_cells, old_y_cells = [a.flatten() for a in np.split(pop.get_cells(), 2, axis = 1)]

    # choose direction using movement surface, if applicable
    if params['movement_surf'] == True:

        mu_distance = params['mu_distance']
        sigma_distance = params['sigma_distance']

        #create a vector of random choices from the approximation vectors for the VonMises mixture dist. for each individual's movement_surface cell
        choices = r.randint(low = 0, high = land.movement_surf_approx_len, size = len(old_x_cells))
        #and use those choices to draw movement directions
        direction = land.movement_surf[old_y_cells, old_x_cells, choices]


        # NOTE: Pretty sure that I don't need to constrain values output for the Gaussian KDE that is approximating the von Mises mixture distribution to 0<=val<=2*pi, because e.g. cos(2*pi + 1) = cos(1), etc...
        # NOTE: indexed out of movement_surf as y then x because becuase the list of lists (like a numpy array structure) is indexed i then j, i.e. vertical, then horizontal



    # else, choose direction using a random walk with a uniform vonmises
    elif params['movement_surf'] == False:

        mu_direction = params['mu_direction']
        kappa_direction = params['kappa_direction']
        mu_distance = params['mu_distance']
        sigma_distance = params['sigma_distance']

        direction = r_vonmises(mu_direction, kappa_direction, size = len(old_x))

    # choose distance
    # NOTE: Instead of lognormal, could use something with long right tail for Levy-flight type movement, same as below
    distance = wald(mu_distance, sigma_distance, size = len(old_x))
    # distance = lognormal(mu_distance, sigma_distance)
    # distance = gamma(mu_distance, sigma_distance)
    #new_pos = (min(max(individual.x + cos(direction) * distance, 0), land.dims[0] - 0.001),
              #min(max(individual.y + sin(direction) * distance, 0), land.dims[1] - 0.001))

    #create the new locations by adding x- and y-dim line segments to their current positions, using trig
    new_x = old_x + cos(direction)*distance
    #and clip the values to be within the landscape dimensions
        #NOTE: subtract a very small value to avoid having the dimension itself set as a coordinate, which rounds down to a cell id one beyond the largest cell id the landscape
    new_x = np.clip(new_x, a_min = 0, a_max = land.dims[1]-0.00001) 
    new_y = old_y + sin(direction)*distance
    new_y = np.clip(new_y, a_min = 0, a_max = land.dims[0]-0.00001)

    #then feed the new locations into each individual's set_pos method
    [ind.set_pos(new_x[n], new_y[n]) for n, ind in enumerate(pop.individs.values())]
    #land.mating_grid.move(old_pos, new_pos, individual)






# Function to generate a simulative Von Mises mixture distribution sampler function

# Returns a lambda function that is a quick and reliable way to simulate draws from a Von Mises mixture distribution:
# 1.) Chooses a direction by neighborhood-weighted probability
# 2.) Makes a random draw from a Von Mises dist centered on the direction, with a kappa value set such that the net effect, when doing this a ton of times for a given neighborhood and then plotting the resulting histogram, gives the visually/heuristically satisfying approximation of a Von Mises mixture distribution


# Runs the Von Mises mixture sampler function (create_von_mises_mix_sampler) across the entire landscape and returns an array-like (list of
# lists) of the resulting lambda-function samplers

def create_movement_surface(land, params, kappa=12, approx_len = 5000):
    queen_dirs = np.array([[-3 * pi / 4, -pi / 2, -pi / 4], [pi, np.NaN, 0], [3 * pi / 4, pi / 2, pi / 4]])

    # grab the correct landscape raster
    rast = land.scapes[params['n_movement_surf_scape']].raster.copy()

    # create embedded raster (so that the edge probabilities are appropriately calculated)
    embedded_rast = np.zeros(shape=[n + 2 for n in rast.shape])
    embedded_rast[1:embedded_rast.shape[0] - 1, 1:embedded_rast.shape[1] - 1] = rast

    # create list of lists (aping an array) for storage of resulting functions
    #movement_surf = [[None for j in range(land.dims[1])] for i in range(land.dims[0])]
    #NOTE: nevermind that, create an actual array and store vectors approximating the functions!
    movement_surf = np.float16(np.zeros((rast.shape[0], rast.shape[1], approx_len)))

    for i in range(rast.shape[0]):
        for j in range(rast.shape[1]):
            neigh = embedded_rast[i:i + 3, j:j + 3].copy()
            #movement_surf[i][j] = create_von_mises_mix_sampler(neigh, queen_dirs, kappa=kappa)
            movement_surf[i, j, :] = create_von_mises_mix_sampler(neigh, queen_dirs, kappa = kappa, approx_len= approx_len)

    return (movement_surf)


def create_von_mises_mix_sampler(neigh, dirs, kappa=12, approx_len = 5000):
    # NOTE: Just visually, heuristically, kappa = 10 seemed like a perfect middle value (kappa ~3 gives too
    # wide of a Von Mises variance and just chooses values around the entire circle regardless of neighborhood
    # probability, whereas kappa ~20 produces noticeable reductions in probability of moving to directions
    # between the 8 queen's neighborhood directions (i.e. doesn't simulate the mixing well enough), and would
    # generate artefactual movement behavior); 12 also seemed to really well in generating probability valleys
    # when tested on neighborhoods that should generate bimodal distributions
    d = list(dirs.ravel())
    n = list(neigh.ravel())
    del d[4]
    del n[4]
    sum_n = float(sum(n))
    n = [i / sum_n for i in n]
    s_vonmises.a = -np.inf
    s_vonmises.b = np.inf
    loc_choices = r.choice(d, approx_len, replace = True, p = n)
    loc_choices = list(C(loc_choices).items())
    approx = np.hstack([s_vonmises.rvs(kappa, loc=loc, scale=1, size = size) for loc, size in loc_choices])
    return(approx)


# function for plotting average unit vectors across the movement surface, for visualization (and for debugging the movement surface functions)
def plot_movement_surf_vectors(land, params, circle=False):
    from numpy import pi, mean, cos, sin, arctan
    from math import sqrt

    # plot movement surface raster
    land.show(scape_num=params['n_movement_surf_scape'])

    # define inner function for plotting a single cell's average unit vector
    def plot_one_cell(i, j):
        # draw sample of angles from the Gaussian KDE representing the von mises mixture distribution (KDE)
        samp = random.sample(list(land.movement_surf[i][j]), k = 100)

        # create lists of the x and y (i.e. cos and sin) components of each angle in the sample
        x_vects = cos(samp)
        y_vects = sin(samp)

        # define x and y plotting coordinates for base of arrow
        # NOTE: they are just equal to the cell's j,i indicies, because I would add 0.5 to plot base of the arrow in the cell center
        # but then would subtract 0.5 to visually reconcile the offset between the plotting axes and the raster
        x, y = j, i

        if circle:  # NOTE: This was just an offhand thought while I was waiting for something else to
            # compute. Doesn't work at all yet, but would be nice to get working. (Would plot circular
            # distributions centered in each cell and scaled to unity, instead of the average vectors I
            # currently have it plotting.)
            x += j
            y += i

            plt.plot(x, y, '.', color='black')


        else:
            # define the dx and dy distances used to the position the arrowhead
            # NOTE: multiply by sqrt(2)/2, to scale to the size of half of the diagonal of a cell
            dx = mean(x_vects) / sqrt(2)
            dy = mean(y_vects) / sqrt(2)
            # NOTE: need to invert the dy value, so that it will plot correctly on the inverted axes (remember that the axes are inverted because the raster is plotted using imshow, which displays raster rows starting with row 0 at the top and working downward)
            # dy *= -1

            # now plot the arrow
            plt.arrow(x, y, dx, dy, alpha=0.75, color='black', head_width=0.24, head_length=0.32)

    # call the internally defined function as a nested list comprehension for all raster cells, which I believe should do its best to vectorize the whole operation
    [[plot_one_cell(i, j) for i in range(len(land.movement_surf))] for j in range(len(land.movement_surf[0]))]






def plot_scaling_random_choice_with_array_size():
    sizes = [1e2, 1e3, 1e4, 1e5, 1e6]
    mean_times = []
    for size in sizes:
        arr = r.random(size = int(size))
        times = []
        for i in range(10):
            start = time.time()
            r.choice(arr)
            stop = time.time()
            diff = stop-start
            times.append(diff)
        mean_times.append(mean(times))
    plt.plot(sizes, mean_times)




