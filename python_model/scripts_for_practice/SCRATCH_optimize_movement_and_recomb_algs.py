#!/usr/bin/python

'''
Just random work I'm doing while poking at and trying optimize the movement and recombination algorithms.

NOTES: In this docstring

CODE: Below


#########
MOVEMENT:
#########

- check how r.choice(n) scales with length of n
- check %lprof for both


N         pop.move()   new_surf_move()  
1000      139ms        4.35ms
10000     1.41s        1.14s



ANOTHER COMPARSION USING %timeit (USING MOVEMENT SURFACE APPROACHES BELOW AND COMPARING TO ORIGINAL):

N               orig                    alt_movement_surf()         alt_alt_movement_surf()
1000            214+-3.87 ms            9.28 ms +- 80.8 us          1.63 ms +- 36.5 us
10000           1.08 s +- 11.2 ms       46.3 ms +- 727 us           7.97 ms +- 51.7 us

AND COMPARING THE main() FN ON A DEMOGRAPHICALLY STABLE POP, USING THE SAME 3 DIFFERENT MOVE ALGS (A
N~1000          1.45 s +- 58.9 ms       1.21 s +- 58.9 ms           1.27 s +- 71.1 ms


Why should it seem to perform so much better on its own, but make less difference in the larger algorithm? I
guess it IS just because movement is only a portion of the total runtime.A

When I run the main() function with the three above algorithms and a pop ~1000 I get:


orig                       alt_movement_surf()         alt_alt_movement_surf()
7.683 s total time         6.006 s total time          5.801 s total time 
movement is 28.3% of time  movement is 1.3% of time    movement is 0.4% of time

##############
RECOMBINATION:
##############

Maybe iteratively convert and store chunks of the total recombination approx array as bits!?! then just return
a random bit from each inter-locus position





- read about and play with int types for vector of binomial-draw approx in numpy, for recombination





'''

import numpy as np
import numpy.random as r
import matplotlib as mpl
import matplotlib.pyplot as plt
import time, os, sys
import copy
from operator import itemgetter

import numpy as np
from numpy import sin, cos
from numpy import pi
import numpy.random as r
from numpy.random import vonmises as r_vonmises
from numpy.random import lognormal, choice, gamma, wald
import matplotlib.pyplot as plt
from scipy.stats import vonmises as s_vonmises
from numba import jit

s_vonmises.a = -np.inf
s_vonmises.b = np.inf




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


            
def alt_move(pop, land, params):

    #old_pos = (individual.x, individual.y)
    ids = pop.individs.keys()
    #print(ids)
    ind_ig = itemgetter(*ids)
    old_x = np.array(ind_ig(pop.get_x_coords()))
    old_y = np.array(ind_ig(pop.get_y_coords()))

    # choose direction using movement surface, if applicable
    if params['movement_surf'] == True:

        cell_strings = landscape.get_cell_strings(old_y, old_x, land.dim_om)
        ig = itemgetter(*cell_strings)


        mu_distance = params['mu_distance']
        sigma_distance = params['sigma_distance']

        # OLD NOTE: this should be a function that allows individuals' movements to be determined by some calculation on a resistance surface raster yet to figure this out; was thinking perhaps of surrounding the surface with a ring of cells with value 0, then calculating a Von Mises mixture distribution (a mixture of all 8 queen's neighbourhood cells) for each cell, then using that in each cell to choose the direction of an individual's movement...

        # NOTE: DEH 01-05-18: For some reason, this line produced a numeric on the office computer, but a numpy array on my laptop. So wrapping it in np.atleast_1d() before indexing [0] should make it work on both machines. To really 'fix' it, might need to upgrade the package on one of the two machines though...

        #replacing list lookup with cell-strings method
        #direction = land.movement_surf[int(individual.y)][int(individual.x)]()
        direction = np.array([r.choice(f) for f in ig(land.alt_movement_surf)])


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

    new_x = old_x + cos(direction)*distance
    new_x = np.clip(new_x, a_min = 0, a_max = land.dims[1]-0.00001)
    new_y = old_y + sin(direction)*distance
    new_y = np.clip(new_y, a_min = 0, a_max = land.dims[0]-0.00001)

    [pop.individs[ind].set_pos(new_x[i], new_y[i]) for i, ind in enumerate(ids)]
    #land.mating_grid.move(old_pos, new_pos, individual)


# Function to generate a simulative Von Mises mixture distribution sampler function

# Returns a lambda function that is a quick and reliable way to simulate draws from a Von Mises mixture distribution:
# 1.) Chooses a direction by neighborhood-weighted probability
# 2.) Makes a random draw from a Von Mises dist centered on the direction, with a kappa value set such that the net effect, when doing this a ton of times for a given neighborhood and then plotting the resulting histogram, gives the visually/heuristically satisfying approximation of a Von Mises mixture distribution

       


def alt_alt_move(pop, land, params):

    #old_pos = (individual.x, individual.y)
    ids = pop.individs.keys()
    #print(ids)
    ind_ig = itemgetter(*ids)
    old_x = np.array(ind_ig(pop.get_x_coords()))
    old_y = np.array(ind_ig(pop.get_y_coords()))
    old_x_cells = np.int16(old_x)
    old_y_cells = np.int16(old_y)

    # choose direction using movement surface, if applicable
    if params['movement_surf'] == True:

        mu_distance = params['mu_distance']
        sigma_distance = params['sigma_distance']

        # OLD NOTE: this should be a function that allows individuals' movements to be determined by some calculation on a resistance surface raster yet to figure this out; was thinking perhaps of surrounding the surface with a ring of cells with value 0, then calculating a Von Mises mixture distribution (a mixture of all 8 queen's neighbourhood cells) for each cell, then using that in each cell to choose the direction of an individual's movement...

        # NOTE: DEH 01-05-18: For some reason, this line produced a numeric on the office computer, but a numpy array on my laptop. So wrapping it in np.atleast_1d() before indexing [0] should make it work on both machines. To really 'fix' it, might need to upgrade the package on one of the two machines though...

        #replacing list lookup with cell-strings method
        #direction = land.movement_surf[int(individual.y)][int(individual.x)]()
        choices = r.randint(low = 0, high = 500, size = len(old_x_cells))
        direction = land.alt_alt_movement_surf[old_y_cells, old_x_cells, choices]


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

    new_x = old_x + cos(direction)*distance
    new_x = np.clip(new_x, a_min = 0, a_max = land.dims[1]-0.00001)
    new_y = old_y + sin(direction)*distance
    new_y = np.clip(new_y, a_min = 0, a_max = land.dims[0]-0.00001)

    [pop.individs[ind].set_pos(new_x[i], new_y[i]) for i, ind in enumerate(ids)]
    #land.mating_grid.move(old_pos, new_pos, individual)






def alt_create_movement_surface(land, params, kappa=12, approx_len = 500):
    # queen_dirs = np.array([[5*pi/4, 3*pi/2, 7*pi/4],[pi, np.NaN, 0],[3*pi/4,pi/2,pi/4]])
    queen_dirs = np.array([[-3 * pi / 4, -pi / 2, -pi / 4], [pi, np.NaN, 0], [3 * pi / 4, pi / 2, pi / 4]])

    # support = np.linspace(s_vonmises.ppf(10e-13, 3, loc = 0), s_vonmises.ppf(1-10e-13, 3, loc = 0), 100000)

    # grab the correct landscape raster
    rast = land.scapes[params['n_movement_surf_scape']].raster.copy()

    # create embedded raster (so that the edge probabilities are appropriately calculated)
    embedded_rast = np.zeros(shape=[n + 2 for n in rast.shape])
    embedded_rast[1:embedded_rast.shape[0] - 1, 1:embedded_rast.shape[1] - 1] = rast

    # create list of lists (aping an array) for storage of resulting functions
    #movement_surf = [[None for j in range(land.dims[1])] for i in range(land.dims[0])]

    gi, gj = np.mgrid[0:land.dims[1], 0:land.dims[0]]
    gi = gi.ravel()
    gj = gj.ravel()
    cells = landscape.get_cell_strings(gi, gj, max([len(str(land.dims[0])), len(str(land.dims[1]))]))
    cells = cells[::-1]

    movement_surf = {}

    for i in range(land.dims[0]):
        for j in range(land.dims[1]):
            neigh = embedded_rast[i:i + 3, j:j + 3].copy()
            movement_surf[cells.pop()] = alt_gen_von_mises_mix_sampler(neigh, queen_dirs, kappa=kappa, approx_len = approx_len)

    return (movement_surf)


def alt_alt_create_movement_surface(land, params, kappa=12, approx_len=500):
    # queen_dirs = np.array([[5*pi/4, 3*pi/2, 7*pi/4],[pi, np.NaN, 0],[3*pi/4,pi/2,pi/4]])
    queen_dirs = np.array([[-3 * pi / 4, -pi / 2, -pi / 4], [pi, np.NaN, 0], [3 * pi / 4, pi / 2, pi / 4]])

    # support = np.linspace(s_vonmises.ppf(10e-13, 3, loc = 0), s_vonmises.ppf(1-10e-13, 3, loc = 0), 100000)

    # grab the correct landscape raster
    rast = land.scapes[params['n_movement_surf_scape']].raster.copy()

    # create embedded raster (so that the edge probabilities are appropriately calculated)
    embedded_rast = np.zeros(shape=[n + 2 for n in rast.shape])
    embedded_rast[1:embedded_rast.shape[0] - 1, 1:embedded_rast.shape[1] - 1] = rast

    # create list of lists (aping an array) for storage of resulting functions
    #movement_surf = [[None for j in range(land.dims[1])] for i in range(land.dims[0])]

    gi, gj = np.mgrid[0:land.dims[1], 0:land.dims[0]]
    gi = gi.ravel()
    gj = gj.ravel()
    cells = landscape.get_cell_strings(gi, gj, max([len(str(land.dims[0])), len(str(land.dims[1]))]))
    cells = cells[::-1]

    movement_surf = np.zeros((land.dims[0], land.dims[1], approx_len))

    for i in range(land.dims[0]):
        for j in range(land.dims[1]):
            neigh = embedded_rast[i:i + 3, j:j + 3].copy()
            movement_surf[i, j, :] = alt_gen_von_mises_mix_sampler(neigh, queen_dirs, kappa=kappa, approx_len = approx_len)

    return (movement_surf)









# Function to generate a simulative Von Mises mixture distribution sampler function

# Returns a lambda function that is a quick and reliable way to simulate draws from a Von Mises mixture distribution:
# 1.) Chooses a direction by neighborhood-weighted probability
# 2.) Makes a random draw from a Von Mises dist centered on the direction, with a kappa value set such that the net effect, when doing this a ton of times for a given neighborhood and then plotting the resulting histogram, gives the visually/heuristically satisfying approximation of a Von Mises mixture distribution

def alt_gen_von_mises_mix_sampler(neigh, dirs, kappa=12, approx_len = 500):
    # NOTE: Just visually, heuristically, kappa = 10 seemed like a perfect middle value (kappa ~3 gives too
    # wide of a Von Mises variance and just chooses values around the entire circle regardless of neighborhood
    # probability, whereas kappa ~20 produces noticeable reductions in probability of moving to directions
    # between the 8 queen's neighborhood directions (i.e. doesn't simulate the mixing well enough), and would
    # generate artefactual movement behavior); 12 also seemed to really well in generating probability valleys
    # when tested on neighborhoods that should generate bimodal distributions
    d = list(dirs.ravel())
    n = list(neigh.ravel())
    # print(d)
    # print(n)
    del d[4]
    del n[4]
    sum_n = float(sum(n))
    n = [i / sum_n for i in n]
    s_vonmises.a = -np.inf
    s_vonmises.b = np.inf
    f = lambda: s_vonmises.rvs(kappa, loc=r.choice(d, 1, p=n), scale=1)
    approx = np.array([f() for _ in range(approx_len)]).ravel()
    f = approx
    return (f)

    # weighted_pdfs = np.array([s_vonmises.pdf(support[i,:], kappa, scale = 1, loc = d[i])*(n[i]/float(sum(n))).ravel() for i in range(len(d))])
    # weighted_pdfs = np.array([s_vonmises.pdf(support, kappa, scale = 1, loc = d[i])*(n[i]/float(sum(n))).ravel() for i in range(len(d))])

    # mix_pdf = weighted_pdfs.sum(axis = 0)

    # return(mix_pdf)






#COULD I USE THIS IN SOME WAY TO SPEED UP DRAWS FROM RANDOM NUMBERS?
    #vonmises approximations could be land.dims[0] x land.dims[1] x ~1e5 arrays that could be subsetted like this perhaps?
def compare_two_approaches_drawing_random_nums():
    test = r.random(size = 100000).reshape((1000,100))
    test2 = [test[i,:] for i in test.shape[0]]
    test[:,r.randint(low = 0, high = test.shape[1], size = test.shape[0])]
    [r.choice(t) for t in test2] #%timeit calls on both reveal that this one is 1.5x the runtime of the previous line


