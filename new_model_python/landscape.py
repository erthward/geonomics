#!/usr/bin/python
#landscape.py

'''
##########################################

Module name:          landscape


Module contains:
                      - definition of the Landscape type
                      - function for creating a random landscape, based on input parameters
                      - associated functions


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           12-28-15
Documentation:        URL


##########################################
'''

from scipy import interpolate
import numpy as np
import numpy.random as r
import matplotlib as mpl



#------------------------------------
# CLASSES ---------------------------
#------------------------------------


class Landscape:
    def __init__(self, dims, raster):
        self.dims = dims
        self.raster = raster
        assert type(self.dims) in [tuple, list], "dims must be expressed on a tuple or a list"
        assert type(self.raster) == np.ndarray, "raster should be a numpy.ndarray"


    #####################
    ### OTHER METHODS ###
    #####################


    def show(self, colorbar = True):
        cmap = 'terrain'
        mpl.pyplot.imshow(self.raster, interpolation = 'nearest', cmap = cmap)
        if colorbar:
            mpl.pyplot.colorbar()




class Landscape_Stack([random_surface(dims, n_rand_pts, interp_method = interp_method[n], num_hab_types = num_hab_types) for n in range(num_scapes)])




