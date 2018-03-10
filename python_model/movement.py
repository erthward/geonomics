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
from numpy import pi
import numpy.random as r
from numpy.random import vonmises as r_vonmises
from numpy.random import lognormal, choice, gamma, wald
import matplotlib.pyplot as plt
from scipy.stats import vonmises as s_vonmises

s_vonmises.a = -np.inf
s_vonmises.b = np.inf




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
       

		#NOTE: DEH 01-05-18: For some reason, this line produced a numeric on the office computer, but a numpy array on my laptop. So wrapping it in np.atleast_1d() before indexing [0] should make it work on both machines. To really 'fix' it, might need to upgrade the package on one of the two machines though...
        direction = np.atleast_1d(land.movement_surf[int(individual.y)][int(individual.x)]())[0]
        #NOTE: Pretty sure that I don't need to constrain values output for the Gaussian KDE that is approximating the von Mises mixture distribution to 0<=val<=2*pi, because e.g. cos(2*pi + 1) = cos(1), etc...
        #NOTE: indexed out of movement_surf as y then x because becuase the list of lists (like a numpy array structure) is indexed i then j, i.e. vertical, then horizontal



    #else, choose direction using a random walk with a uniform vonmises
    elif params['movement_surf'] == False:

        mu_direction = params['mu_direction']
        kappa_direction = params['kappa_direction']
        mu_distance = params['mu_distance']
        sigma_distance = params['sigma_distance']


        direction = r_vonmises(mu_direction, kappa_direction)

   
    #choose distance
    #NOTE: Instead of lognormal, could use something with long right tail for Levy-flight type movement, same as below
    distance = wald(mu_distance, sigma_distance)
    #distance = lognormal(mu_distance, sigma_distance) 
    #distance = gamma(mu_distance, sigma_distance) 


    #determine new coordinates
    individual.x = min(max(individual.x + cos(direction)*distance, 0), land.dims[0]-0.001)
    individual.y = min(max(individual.y + sin(direction)*distance, 0), land.dims[1]-0.001)





#Function to generate a simulative Von Mises mixture distribution sampler function
    
    #Returns a lambda function that is a quick and reliable way to simulate draws from a Von Mises mixture distribution:
    #1.) Chooses a direction by neighborhood-weighted probability
    #2.) Makes a random draw from a Von Mises dist centered on the direction, with a kappa value set such that the net effect, when doing this a ton of times for a given neighborhood and then plotting the resulting histogram, gives the visually/heuristically satisfying approximation of a Von Mises mixture distribution

def gen_von_mises_mix_sampler(neigh, dirs, kappa=12): 
    #NOTE: Just visually, heuristically, kappa = 10 seemed like a perfect middle value (kappa ~3 gives too
    #wide of a Von Mises variance and just chooses values around the entire circle regardless of neighborhood
    #probability, whereas kappa ~20 produces noticeable reductions in probability of moving to directions
    #between the 8 queen's neighborhood directions (i.e. doesn't simulate the mixing well enough), and would
    #generate artefactual movement behavior); 12 also seemed to really well in generating probability valleys
    #when tested on neighborhoods that should generate bimodal distributions
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
    f = lambda: s_vonmises.rvs(kappa, loc = r.choice(d, 1, p = n), scale = 1)
    return(f)
    

    #weighted_pdfs = np.array([s_vonmises.pdf(support[i,:], kappa, scale = 1, loc = d[i])*(n[i]/float(sum(n))).ravel() for i in range(len(d))])
    #weighted_pdfs = np.array([s_vonmises.pdf(support, kappa, scale = 1, loc = d[i])*(n[i]/float(sum(n))).ravel() for i in range(len(d))])


    #mix_pdf = weighted_pdfs.sum(axis = 0)

    #return(mix_pdf)




#Runs the Von Mises mixture sampler function (gen_von_mises_mix_sampler) across the entire landscape and returns an array-like (list of
#lists) of the resulting lambda-function samplers

def create_movement_surface(land, params, kappa = 12):

    #queen_dirs = np.array([[5*pi/4, 3*pi/2, 7*pi/4],[pi, np.NaN, 0],[3*pi/4,pi/2,pi/4]])
    queen_dirs = np.array([[-3*pi/4,-pi/2,-pi/4],[pi, np.NaN, 0],[3*pi/4,pi/2,pi/4]])

    #support = np.linspace(s_vonmises.ppf(10e-13, 3, loc = 0), s_vonmises.ppf(1-10e-13, 3, loc = 0), 100000)


    #grab the correct landscape raster
    rast = land.scapes[params['n_movement_surf_scape']].raster.copy()

    #create embedded raster (so that the edge probabilities are appropriately calculated)
    embedded_rast = np.zeros(shape = [n+2 for n in rast.shape])
    embedded_rast[1:embedded_rast.shape[0]-1, 1:embedded_rast.shape[1]-1] = rast
    

    #create list of lists (aping an array) for storage of resulting functions
    movement_surf = [[None for j in range(land.dims[1])] for i in range(land.dims[0])]

    for i in range(rast.shape[0]):
        for j in range(rast.shape[1]):
            neigh = embedded_rast[i:i+3, j:j+3].copy()
            movement_surf[i][j] = gen_von_mises_mix_sampler(neigh, queen_dirs, kappa = kappa)


    return(movement_surf)







#function for plotting average unit vectors across the movement surface, for visualization (and for debugging the movement surface functions)
def plot_movement_surf_vectors(land, params, circle = False):
    from numpy import pi, mean, sqrt, cos, sin, arctan

    #plot movement surface raster
    land.show(scape_num = params['n_movement_surf_scape'])

    #define inner function for plotting a single cell's average unit vector
    def plot_one_cell(i,j):
        #draw sample of angles from the Gaussian KDE representing the von mises mixture distribution (KDE)
        #jason: changed from 1000 to 100 
        samp =  np.array([land.movement_surf[i][j]() for n in range(100)])

        #create lists of the x and y (i.e. cos and sin) components of each angle in the sample
        x_vects = cos(samp)
        y_vects = sin(samp)
           


        #define x and y plotting coordinates for base of arrow
            #NOTE: they are just equal to the cell's j,i indicies, because I would add 0.5 to plot base of the arrow in the cell center
            #but then would subtract 0.5 to visually reconcile the offset between the plotting axes and the raster
        x,y = j,i

        if circle: #NOTE: This was just an offhand thought while I was waiting for something else to
            #compute. Doesn't work at all yet, but would be nice to get working. (Would plot circular
            #distributions centered in each cell and scaled to unity, instead of the average vectors I
            #currently have it plotting.)
            x += j
            y += i

            plt.plot(x,y, '.', color = 'black')


        else:
            #define the dx and dy distances used to the position the arrowhead
            #NOTE: multiply by sqrt(2)/2, to scale to the size of half of the diagonal of a cell
            dx = mean(x_vects)/sqrt(2)
            dy = mean(y_vects)/sqrt(2)       
            #NOTE: need to invert the dy value, so that it will plot correctly on the inverted axes (remember that the axes are inverted because the raster is plotted using imshow, which displays raster rows starting with row 0 at the top and working downward)
            #dy *= -1
       

            #now plot the arrow
            plt.arrow(x, y, dx, dy, alpha = 0.75, color = 'black', head_width = 0.24, head_length = 0.32)


    #call the internally defined function as a nested list comprehension for all raster cells, which I believe should do its best to vectorize the whole operation
    [[plot_one_cell(i,j) for i in range(len(land.movement_surf))] for j in range(len(land.movement_surf[0]))]







