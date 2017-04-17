#!/usr/bin/python
#NEW_pop_dynamics.py

import numpy as np
import landscape
from scipy import interpolate

'''A totally new take on population dynamics.'''



def calc_pop_dens(land, pop, window_width = None, normalize_by = 'census', max_1 = False):

    '''
    Calculate an interpolated raster of local population density, using a window size of window_width.
    Valid values for normalize_by currently include 'census' and 'none'. If normalize_by = 'census', max_1 =
    True will cause the output density raster to vary between 0 and 1, rather than between 0 and the current
    max normalized density value. Window width will default to 1/10 of the larger raster dimension.
    '''

    #window length default to 1/10 the maximum landscape dimension
    if window_width == None:
        window_width = max(land.dims)*0.1

    #create raster for storing local density vals
    local_dens_rast = np.array([[0]*land.dims[0]]*land.dims[1])

    #shorthand
    dims = land.dims

    #get a list of pop's coord-tuples
    c = pop.get_coords().values() 

    #make window_width a float, to avoid Py2 integer-division issues
    window_width = float(window_width)
    
    #create meshgrid using window_width/2 as step size
    grid_j, grid_i = np.mgrid[1:dims[0]:complex("%ij" % (dims[0]/(window_width/2))), 1:dims[1]:complex("%ij" % (dims[1]/(window_width/2)))]

    #flatten the arrays, so that I can run over them in a single for loop
    gj = grid_j.ravel()
    gi = grid_i.ravel()

    #make lists of tuples, of same length as gj, containing the window ll and ur coords
    window_ll = [(max(gj[n]-window_width/2, 0), max(gi[n]-window_width/2, 0)) for n in range(len(gj))]   #constrain min window vals to 0
    window_ur = [(min(gj[n]+window_width/2, land.dims[0]), min(gi[n]+window_width/2, land.dims[1])) for n in range(len(gj))] #constrain max window vals to each respective land dimension
    assert len(window_ll) == len(gj)
    assert len(window_ur) == len(gj)

    #make a list of the sizes of each window
    window_size = [(window_ur[n][0] - window_ll[n][0]) * (window_ur[n][1] - window_ll[n][1]) for n in range(len(gj))]#calculate size of this window (not always the same because of edge effects
    assert len(window_size) == len(gj)
    
    #make a list of the counts of all individs within each window
    window_ct = [len([ind for ind in range(len(c)) if (c[ind][0]>window_ll[n][0] and c[ind][0]<=window_ur[n][0]) and (c[ind][1]>window_ll[n][1] and c[ind][1]<=window_ur[n][1])]) for n in range(len(gj))] 
    assert len(window_ct) == len(gj)

    #divide window counts by window sizes
    window_dens = [window_ct[n]/window_size[n] for n in range(len(window_ct))] #divide by window size
    assert len(window_dens) == len(gj)

    #if normalize_by == census, then divide each density by total pop census size
    if normalize_by == 'census':
        N = pop.census()
        window_dens = [dens/N for dens in window_dens]
    elif normalize_by == 'none':
        pass

    else:  #POTENTIALLY ADD OTHER OPTIONS HERE, i.e. to normalize by starting size?
        pass 

    #interpolate resulting density vals to a grid equal in size to the landscape
    new_gj, new_gi = np.mgrid[1:dims[0]:complex("%ij" % (dims[0])), 1:dims[1]:complex("%ij" % (dims[1]))]
    dens = interpolate.griddata(zip(list(gj), list(gi)), window_dens, (new_gj, new_gi), method = 'cubic')



    if normalize_by <> 'none':

        #if max_1 == True, set max_val to dens.max(), such that the density raster output will be normalized to
        #its own max, and thus vary between 0 and 1; else set to 1, and the output raster will vary between 0 and the current max value
        if max_1 == True:
            max_val = dens.max()
        elif max_1 == False:
            max_val = 1

        #Use max_val to normalize the density raster to either 0 to its current max val or
        #0 to 1, to make sure the interpolation didn't generate any values slightly outside this range
        norm_factor = max_val - dens.min()
        dens = (dens - dens.min())/norm_factor

    return(landscape.Landscape(dims, dens))



def logistic_growth(r, N, K):
    K = float(K)
    dN_dt = r*(1-(N/K))*N
    return(dN_dt)






def burn_in(pop, land, params):

    #NOTE: Still not totally clear to me how to get this to work, but I want to use some iterative process to
    #get the model to empirically decide its own stable spatial population distribution based on its given
    #parameterization and start conditions, and to break the burn-in only once it has reached that state. I am
    #thinking that to do this I can use the previous timestep's population density as a local carrying capacity
    #for density-dependent births/deaths, and that by iterating the opposing forces of movement
    #to/through high-permeability areas and attrition in those areas because of density-dependence should reach
    #a dynamic equilibrium, which the model can determine through assessment of a metric or metrics




    #start the mean_dens array as just the starting pop density
    mean_dens = calc_pop_dens(land, pop)


    #a while loop conditioned on burn == True, which will only switch to False once the stasis/stationarity criteria are determined to be met
    #while burn == True:
    for t in burn_T:

        #calculate pop density array, which will function as the carrying capacity (K) in each step
        K = calc_pop_dens(land, pop)
        
        #calc local logistic growth rates

        
        #enact births and deaths, on basis of local logistic growth rates
            #NOTE: Probably should not be using selection here, as I want this to determine the initial
            #dynamic-equilibrium distribution on the basis of movement and density-dependence only


        #recalculate mean_dens (using the formula I worked out on the board for updating a mean)
        mean_dens = mean_dens + (1/(t+1))*(K-mean_dens)

        #assess stationarity of pop density
        
        #if determined stationary, based on chosen metric(s), 
            #set burn = False

        #else
            #repeat loop



        




