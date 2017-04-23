#!/usr/bin/python


'''Trying to redo the VonMises mixture thing and make it much smoother and faster.'''


import numpy as np
import numpy.random as r
from numpy import pi
from scipy.stats import vonmises

vonmises.a = -np.inf
vonmises.b = np.inf

def gen_von_mises_mix_sampler(neigh, dirs, kappa=12): 
    #Returns a lambda function that is a quick and reliable way to simulate draws from a Von Mises mixture distribution:
    #1.) Choose a direction by neighborhood-weighted probability
    #2.) Make a random draw from a Von Mises dist centered on the direction, with a kappa value set such that the net effect, when doing this a ton of times for a given neighborhood and then plotting the resulting histogram, gives the visually/heuristically satisfying approximation of a Von Mises mixture distribution

    #NOTE: Just visually, heuristically, kappa = 10 seemed like a perfect middle value (kappa ~3 gives too
    #wide of a Von Mises variance and just chooses values around the entire circle regardless of neighborhood
    #probability, whereas kappa ~20 produces noticeable reductions in probability of moving to directions
    #between the 8 queen's neighborhood directions (i.e. doesn't simulate the mixing well enough), and would
    #generate artefactual movement behavior); 12 also seemed to really well in generating probability valleys
    #when tested on neighborhoods that should generate bimodal distributions
    d = list(dirs.ravel())
    n = list(neigh.ravel())
    del d[4]
    del n[4]
    n = [i/float(sum(n)) for i in n]
    vonmises.a = -np.inf
    vonmises.b = np.inf
    f = lambda: vonmises.rvs(kappa, loc = r.choice(d, 1, p = n), scale = 1)[0]
    return(f)
    

    #weighted_pdfs = np.array([vonmises.pdf(support[i,:], kappa, scale = 1, loc = d[i])*(n[i]/float(sum(n))).ravel() for i in range(len(d))])
    #weighted_pdfs = np.array([vonmises.pdf(support, kappa, scale = 1, loc = d[i])*(n[i]/float(sum(n))).ravel() for i in range(len(d))])


    #mix_pdf = weighted_pdfs.sum(axis = 0)

    #return(mix_pdf)






def create_movement_surface(land, params, kappa = 12):

    #queen_dirs = np.array([[5*pi/4, 3*pi/2, 7*pi/4],[pi, np.NaN, 0],[3*pi/4,pi/2,pi/4]])
    queen_dirs = np.array([[-3*pi/4,-pi/2,-pi/4],[pi, np.NaN, 0],[3*pi/4,pi/2,pi/4]])

    #support = np.linspace(vonmises.ppf(10e-13, 3, loc = 0), vonmises.ppf(1-10e-13, 3, loc = 0), 100000)


    #grab the correct landscape raster
    rast = land.scapes[params['movement_surf_scape_num']].raster.copy()

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



