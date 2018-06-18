#scratch script used to develop the von mises mixture distribution array approach for implementing resistance/movement surfaces

import numpy as np
import numpy.random as r
from sklearn.neighbors.kde import KernelDensity
from scipy.stats import vonmises, vonmises_line




#########
#NOTE: THIS FUNCTION IS UNNECESSARY
#function to get probabilities associated with each queen's neighborhood cell
def probs(neigh):
    probs = r.rand(3,3)*0
    for i in range(3):
        for j in range(3):
            probs[i,j] = neigh[i,j]/(sum(neigh))# - neigh[1,1])  #NOTE: no need to subtract neigh[1,1] bc == 0, from NOTE above
    probs[1,1] = 0
    return probs

#########




#create vector of directions reflecting relative probabilities of each queen's-neighborhood cell, for use in KDE vonmises mixture approximation
def create_representative_dirs_vector(dirs, probs):
    dirs_vector = []
    for i in range(3):
        for j in range(3):
            dirs_vector.extend([dirs[i,j]]*int(round(probs[i,j]*100)))
    return dirs_vector





#roundabout but for now best approach to a function for creating a vonmises mixture distribution 
    #(a la http://stackoverflow.com/questions/28839246/scipy-gaussian-kde-and-circular-data)
def create_approx_vonmises_mixture_dist_KDE(queen_dirs, neigh, kappa, bandwidth, plot = None, superfluous_output = None):
    from scipy.interpolate import interp1d


    #create vector of dirs from queen's neighborhood of cells
    data = create_representative_dirs_vector(queen_dirs, neigh)


    #set limits for vonmises
    vonmises.a = 0
    vonmises.b = 2*np.pi
    x_data = np.linspace(0,2*np.pi, 181) #radian intervals equiv to 2 degrees

    kernels = []

    for d in data:
        kernel = vonmises(kappa, loc=d)
        kernel = kernel.pdf(x_data)
        kernels.append(kernel)

        if plot == True:
            kernel /= kernel.max()
            kernel *= .2
            plt.plot(x_data, kernel, "grey", alpha = 0.5)

    vonmises_kde = np.sum(kernels,axis = 0)
    vonmises_kde = vonmises_kde / np.trapz(vonmises_kde, x = x_data)
    f = interp1d(x_data, vonmises_kde)
    
    #now build large sample drawn from the KDE (i.e. vonmises mixture dist), 
        #in lieu of a pseudorandom number generator based on the KDE, which I haven't yet figured out how to generate...
        #a la http://stats.stackexchange.com/questions/82797/how-to-draw-random-samples-from-a-non-parametric-estimated-distribution
    samp = [vonmises_line.rvs(kappa, loc = r.choice(data, replace = True)) for i in range(7000)]
    #constrain all values to 0 <= val <= 2*pi
    samp = [s - 2*pi if s > 2*pi else s for s in samp]
    samp = [s + 2*pi if s < 0 else s for s in samp]


    

    #TODO: now figure out best computation/memory trade-off for storing this "pseudo-pseudorandom number generator" for each cell in the landscape, to be used in resistance surface-like movement
    #NOTE: FOR NOW AT LEAST, I CAN JUST ESTIMATE A GAUSSIAN KDE BASED ON samp, CREATED ABOVE, THEN SAVE THAT KDE FOR EACH CELL, TO BE LATER USED TO DRAW A DIRECTION
        #NOTE:  I could even constrain the value drawn to be 0 <= val <= 2*pi, although I don't even think that's
                #necessary b/c vals above or below will just drift 'around the circle', so to speak
    kde = KernelDensity(kernel = 'gaussian', bandwidth = bandwidth)
    fit = kde.fit(samp)

    if plot == True:
        plt.plot(x_data, vonmises_kde, c = 'plum', linewidth = 3)
        hist(samp, normed = True, bins = 100, alpha = 0.5)
        hist(fit.sample(1)[0], normed = True, bins = 100, alpha = 0.5)


    if superfluous_output == True:         #if the superfluous tag is True, will vomit all these contents, for exploring/debugging the components
        return f, data, kernels, fit

    else:
        return fit




#function for creating a raster of vonmises mixture approximation KDEs, for use in resistance-surface or habitat-quality movement
def create_vonmises_KDE_array(land, params, plot = None):

    #create array of radian directions corresponding to queen's neighborhood (0rad to right, increasing counter-clockwise) 
    queen_dirs = np.array([[3*pi/4, pi/2, pi/4],[pi, np.NaN, 0],[5*pi/4,3*pi/2,7*pi/4]])

    #copy raster, for use in creating KDE array
    rast = land.scapes[params['movement_scape']].raster.copy()

    #create raster of dimensions 2 greater than rast in each direction, then pad the margin that I added with zeros
    embedded_rast = np.zeros(shape = [n+2 for n in rast.shape])#*np.NaN  
    embedded_rast[1:embedded_rast.shape[0]-1, 1:embedded_rast.shape[1]-1] = rast
    
    #create list of lists, i.e. essentially an array, for storing KDEs
    kde_array= [[None]*rast.shape[1] for i in range(rast.shape[0])]



    #now loop over all cells in rast, but as embedded within nan_rast, and create and store a vonmises mixture approximation KDE for each
    for i in range(rast.shape[0]):
        for j in range(rast.shape[1]):
            neigh = embedded_rast[i:i+3, j:j+3].copy()
            neigh[1,1] = 0
            fit = create_approx_vonmises_mixture_dist_KDE(queen_dirs, neigh, params['movement_scape_vonmises_kappa'], params['movement_scape_gauss_KDE_bandwidth'])
            kde_array[i][j] = fit


    return kde_array









################################################################################################################
#EXAMPLE CODE
#array of radian directions corresponding to queen's neighborhood (0rad to right, increasing counter-clockwise)
L = np.array([[3*pi/4, pi/2, pi/4],[pi, np.NaN, 0],[5*pi/4,3*pi/2,7*pi/4]])

#array of landscape values for queen's neighborhood
l = np.array([[0.5,1,2],[1,0,3],[2,4,5]])/7 #NOTE: l[1,1] == 0


#########
#NOTE: THIS FUNCTION IS REALLY UNNECESSARY
l_probs = probs(l)
#########


kappa = 2
plot = True
superfluous_output = False
fit = create_approx_vonmises_mixture_dist_KDE(L, l, kappa, plot, superfluous_output)
rand_dir = r.choice(fit.sample(1)[0])
#kde_array = create_vonmises_mixture_approx_KDE_array(land, params) #NOTE: this line won't work unless the main_FOR_CMD_LINE.py script is execute, so that a land object and a params dict exist







##############################################################################################################
#NOTE: 

    #SO I THINK THE ABOVE WOULD WORK!
    #BASIC GIST:

    #0.) include params determining what type of movement functions to use, and if this
    #resistance/habitat-quality movement functionality is used, what layer to use to determine
    #(AND IF GENOTYPE AT ANY LOCI SHOULD DETERMINE THE KDE ARRAY THAT APPLIES TO AN INDIVIDUAL??)

    #1.) create KDE of samples drawn from vonmises KDE for each cell (ideally use not just queen's
    #neighborhood but instead a 24-neighbor neighborhood (i.e. the 5x5 neighborhood minus the cell itself), or
    #perhaps even an arbitrarily large neighborhood as determined from comparison between movement distance
    #distribution and cell size?)

    #2.) save all of the KDEs as an array 'underneath' the habitat quality or resistance surface array

    #3.) make quick calls to the KDEs in order to draw random directions for the movement module

