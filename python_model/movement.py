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
from numpy.random import vonmises, lognormal, choice, gamma, wald
import matplotlib.pyplot as plt


#----------------------------------
# TYPES ---------------------------
#----------------------------------






#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------



def move(individual, land, params):
    

    from numpy import sin, cos

    #choose direction using movement surface, if applicable
    if params['movement_surf'] == True:

        mu_distance = params['mu_distance']
        sigma_distance = params['sigma_distance']


        #OLD NOTE: this should be a function that allows individuals' movements to be determined by some calculation on a resistance surface raster yet to figure this out; was thinking perhaps of surrounding the surface with a ring of cells with value 0, then calculating a Von Mises mixture distribution (a mixture of all 8 queen's neighbourhood cells) for each cell, then using that in each cell to choose the direction of an individual's movement...
        
        direction = choice(land.movement_surf[int(np.floor(individual.y))][int(np.floor(individual.x))].sample(1)[0])  
        #NOTE: Pretty sure that I don't need to constrain values output for the Gaussian KDE that is approximating the von Mises mixture distribution to 0<=val<=2*pi, because e.g. cos(2*pi + 1) = cos(1), etc...
        #NOTE: indexed out of movement_surf as y then x because becuase the list of lists (like a numpy array structure) is indexed i then j, i.e. vertical, then horizontal



    #else, choose direction using a random walk with a uniform vonmises
    elif params['movement_surf'] == False:

        mu_direction = params['mu_direction']
        kappa_direction = params['kappa_direction']
        mu_distance = params['mu_distance']
        sigma_distance = params['sigma_distance']


        direction = vonmises(mu_direction, kappa_direction)

   
    #choose distance
    #NOTE: Instead of lognormal, could use something with long right tail for Levy-flight type movement, same as below
    distance = wald(mu_distance, sigma_distance)
    #distance = lognormal(mu_distance, sigma_distance) 
    #distance = gamma(mu_distance, sigma_distance) 


    #determine new coordinates
    new_x = min(max(individual.x + cos(direction)*distance, 0), land.dims[0]-0.001)
    new_y = min(max(individual.y + sin(direction)*distance, 0), land.dims[1]-0.001)

    
    #return new coordinates
    return new_x, new_y




####Set of functions for using the movement surface layer to create an array of Gaussian KDEs approximating
    #von Mises mixture distributions, for implementing non-random probabilistic movement across a resistance-type surface

    #NOTE: BASIC GIST OF THIS APPROACH:

    #0.) include params determining what type of movement functions to use, and if this
    #resistance/habitat-quality movement functionality is used, what layer to use to determine
    #(AND IF GENOTYPE AT ANY LOCI SHOULD DETERMINE THE KDE ARRAY THAT APPLIES TO AN INDIVIDUAL??)

    #1.) create KDE of samples drawn from vonmises KDE for each cell (ideally use not just queen's
    #neighborhood but instead a 24-neighbor neighborhood (i.e. the 5x5 neighborhood minus the cell itself), or
    #perhaps even an arbitrarily large neighborhood as determined from comparison between movement distance
    #distribution and cell size?)

    #2.) save all of the KDEs as an array 'underneath' the habitat quality or resistance surface array

    #3.) make quick calls to the KDEs in order to draw random directions for the movement module





#create vector of directions reflecting relative probabilities of each queen's-neighborhood cell, for use in KDE vonmises mixture approximation
def create_representative_dirs_vector(dirs, probs):
    #print probs
    dirs_vector = []
    for i in range(3):
        for j in range(3):
            dirs_vector.extend([dirs[i,j]]*int(round(probs[i,j]*1000)))
    return dirs_vector





#roundabout but for now best approach to a function for creating a vonmises mixture distribution 
    #(a la http://stackoverflow.com/questions/28839246/scipy-gaussian-kde-and-circular-data)
    #NOTE: I'm leaving the 'plot' and 'verbose' flags in here for now, in case useful for exploration/debugging

def create_approx_vonmises_mixture_dist_KDE(queen_dirs, neigh, kappa, bandwidth, plot = None, verbose = None):
    from numpy import pi
    from scipy.interpolate import interp1d
    from scipy.stats import vonmises, vonmises_line
    from sklearn.neighbors.kde import KernelDensity


    #create vector of dirs from queen's neighborhood of cells
    data = create_representative_dirs_vector(queen_dirs, neigh)
    from collections import Counter as C
    #print len(data)
    #print C(data)


    #set limits for vonmises
    vonmises.a = 0
    vonmises.b = 2*pi
    x_data = np.linspace(0,2*pi, 181) #radian intervals equiv to 2 degrees

    kernels = []

    for d in data:  #TODO:THERE'S LIKELY A MUCH BETTER/QUICKER WAY TO DO THIS THAN BY SUMMING SO MANY PDFs!!!  THIS IS LIKELY THE MAJOR SOURCE OF SLOWNESS
        kernel = vonmises(kappa, loc=d)
        kernel = kernel.pdf(x_data)
        kernels.append(kernel)

        if plot == True:
            kernel /= kernel.max()
            kernel *= .2
            plt.plot(x_data, kernel, "grey", alpha = 0.5)

    vonmises_kde = np.sum(kernels,axis = 0)
    #print vonmises_kde
    vonmises_kde = vonmises_kde / np.trapz(vonmises_kde, x = x_data)
    f = interp1d(x_data, vonmises_kde)
    
    #now build large sample drawn from the KDE (i.e. vonmises mixture dist), 
        #in lieu of a pseudorandom number generator based on the KDE, which I haven't yet figured out how to generate...
        #a la http://stats.stackexchange.com/questions/82797/how-to-draw-random-samples-from-a-non-parametric-estimated-distribution
    samp = [vonmises_line.rvs(kappa, loc = choice(data, replace = True)) for i in range(7000)]
    #constrain all values to 0 <= val <= 2*pi
    samp = [s - 2*pi if s > 2*pi else s for s in samp]
    samp = [s + 2*pi if s < 0 else s for s in samp]


    

    #TODO: now figure out best computation/memory trade-off for storing this "pseudo-pseudorandom number generator" for each cell in the landscape, to be used in resistance surface-like movement
    #NOTE: FOR NOW AT LEAST, I CAN JUST ESTIMATE A GAUSSIAN KDE BASED ON samp, CREATED ABOVE, THEN SAVE THAT KDE FOR EACH CELL, TO BE LATER USED TO DRAW A DIRECTION
    kde = KernelDensity(kernel = 'gaussian', bandwidth = bandwidth)
    fit = kde.fit(samp)

    if plot == True:
        plt.plot(x_data, vonmises_kde, c = 'plum', linewidth = 3)
        hist(samp, normed = True, bins = 100, alpha = 0.5)
        hist(fit.sample(1)[0], normed = True, bins = 100, alpha = 0.5)


    if verbose == True:         #if the verbose tag is True, will vomit all these contents, for exploring/debugging the components
        return f, data, kernels, fit

    else:
        return fit






#function for creating a raster of vonmises mixture approximation KDEs, for use in resistance-surface or habitat-quality movement
def create_vonmises_KDE_array(land, params, plot = None):
    from numpy import pi

    #create array of radian directions corresponding to queen's neighborhood (0rad to right, increasing counter-clockwise) 
        #NOTE: 04-07-17: I BELIEVE I MAY HAVE FOUND THE ROOT OF THE MOVEMENT ISSUE: Radians actually increase
        #clockwise, not counterclockwise, so the way it was set up it was pushing corner individuals into the corners!
    queen_dirs = np.array([[5*pi/4, 3*pi/2, 7*pi/4],[pi, np.NaN, 0],[3*pi/4,pi/2,pi/4]])

    #copy raster, for use in creating KDE array
    rast = land.scapes[params['movement_surf_scape_num']].raster.copy()
    rast[rast == 0] = 10e-4

    #create raster of dimensions 2 greater than rast in each direction, then pad the margin that I added with zeros
    embedded_rast = np.zeros(shape = [n+2 for n in rast.shape])#*np.NaN  
    embedded_rast[1:embedded_rast.shape[0]-1, 1:embedded_rast.shape[1]-1] = rast
    
    #create list of lists, i.e. essentially an array, for storing KDEs
    kde_array= [[None]*rast.shape[1] for i in range(rast.shape[0])]



    #now loop over all cells in rast, but as embedded within nan_rast, and create and store a vonmises mixture approximation KDE for each
    for i in range(rast.shape[0]):
        for j in range(rast.shape[1]):
            #print '\n\n###############\n\t%i, %i\n\n' %(i,j)
            neigh = embedded_rast[i:i+3, j:j+3].copy()
            neigh[1,1] = 0
            fit = create_approx_vonmises_mixture_dist_KDE(queen_dirs, neigh, params['movement_surf_vonmises_kappa'], params['movement_surf_gauss_KDE_bandwidth'])
            kde_array[i][j] = fit


    return kde_array





#function for plotting average unit vectors across the movement surface, for visualization (and for debugging the movement surface functions)
def plot_movement_surf_vectors(land, params):

    from numpy import pi, mean, sqrt, cos, sin, arctan

    #plot movement surface raster
    land.show(scape_num = params['movement_surf_scape_num'])

    #define inner function for plotting a single cell's average unit vector
    def plot_one_cell(i,j):
        #draw sample of angles from the Gaussian KDE representing the von mises mixture distribution (KDE)
        samp = land.movement_surf[i][j].sample()[0]

        #create lists of the x and y (i.e. cos and sin) components of each angle in the sample
        x_vects = [cos(d) for d in samp]
        y_vects = [sin(d) for d in samp]
           


        #define x and y plotting coordinates for base of arrow
            #NOTE: they are just equal to the cell's j,i indicies, because I would add 0.5 to plot base of the arrow in the cell center
            #but then would subtract 0.5 to visually reconcile the offset between the plotting axes and the raster
        x,y = j,i


        #define the dx and dy distances used to the position the arrowhead
            #NOTE: multiply by sqrt(2)/2, to scale to the size of half of the diagonal of a cell
        dx = mean(x_vects)*sqrt(2)/2
        dy = mean(y_vects)*sqrt(2)/2        
        #NOTE: need to invert the dy value, so that it will plot correctly on the inverted axes (remember that the axes are inverted because the raster is plotted using imshow, which displays raster rows starting with row 0 at the top and working downward)
        #dy *= -1
       

        #now plot the arrow
        plt.arrow(x, y, dx, dy, alpha = 0.75, color = 'black', head_width = 0.12, head_length = 0.16)


    #call the internally defined function as a nested list comprehension for all raster cells, which I believe should do its best to vectorize the whole operation
    [[plot_one_cell(i,j) for i in range(len(land.movement_surf))] for j in range(len(land.movement_surf[0]))]


