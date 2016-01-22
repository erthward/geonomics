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


    def show(self):
        cmap = 'terrain'
        mpl.pyplot.imshow(self.raster, interpolation = 'nearest', cmap = cmap)
        mpl.pyplot.colorbar()




class Landscape_stack:
    def __init__(self, raster_list):
        self.scapes = dict(zip(range(len(raster_list)), raster_list))  
        assert False not in [raster.__class__.__name__ == 'Landscape' for raster in raster_list], 'All landscapes supplied in raster_list must be of type landscape.Landscape.'
        self.num_rasters = len(raster_list)
        assert len(set([land.dims for land in self.scapes.values()])) == 1, 'Dimensions of all landscapes must be even.'
        self.dims = self.scapes.values()[0].dims


    def show(self):
        cmaps = ['terrain', 'bone']
        alphas = [1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        for n, scape in self.scapes.items():
            mpl.pyplot.imshow(scape.raster, interpolation = 'nearest', alpha = alphas[n], cmap = cmaps[n])
            mpl.pyplot.colorbar()












#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------


def random_surface(dims, n_rand_pts, interp_method = "cubic", num_hab_types = 2):  #requires landscape to be square, such that dim = domain = range
    #NOTE: can use "nearest" interpolation to create random patches of habitat (by integer values); can change num_hab_types to > 2 to create a randomly multi-classed landscape
    #NOTE: I guess this could be used for rectangular landscapes, if the square raster is generated using the larger of the two dimensions, and then the resulting array is subsetted to the landscape's dimensions
    #NOTE: This seems to generate decent, believable random rasters! 
        # n_rand_pts/dim ratio:
            # ~0.01-0.05 --> broad, simple gradients
            # ~0.05-0.10 --> slightly more complex, occasionally landscape of low regions divided by a major high region (or vice versa)
            # ~0.10-0.50 --> landscape can be broken up into numerous fragmented patches (though sometimes still relatively homogeneous, with one or two small, extremely different patches
            # ~0.50-1.00 --> more highly fragmented version of above
    max_dim = max(dims)
    if interp_method == 'nearest':
        vals = r.rand(n_rand_pts) * (num_hab_types-1)
    else:
        vals = r.rand(n_rand_pts)
    points = r.normal(max_dim/2, max_dim*2,[n_rand_pts,2]) #selects seed points from well outside the eventaul landscape, to ensure interpolation across area of interest
    grid_x, grid_y = np.mgrid[1:max_dim:complex("%ij" % max_dim), 1:max_dim:complex("%ij" % max_dim)]
    I = interpolate.griddata(points, vals, (grid_x, grid_y), method = interp_method)
    if interp_method == 'nearest':  #i.e., if being used to generate random habitat patches...
        I = I.round().astype(int)
    if interp_method == 'cubic':  #transform to constrain all values to 0 <= val <= 1
        I = I + np.abs(I.min())+(0.01*r.rand()) #NOTE: adding a small jitter to keep values from reaching == 0 or == 1, as would likely be the case with linear interpolation
        I = I/(I.max()+(0.01*r.rand()))
    if dims[0] <> dims[1]:
        pass #NOTE: figure out how to get it to use the dims tuple to subset an approriate size if dims not equal
    return Landscape(dims,I)




def build_scape_stack(params, num_hab_types = 2):

    #NOTE: If a multi-class (rather than binary) block-habitat raster would be of interest, would need to make num_hab_types customizable)


    #grab necessary parameters from the params dict

    if params['num_scapes'] == None:
        num_scapes = 1
    else:
        num_scapes = params['num_scapes']


    if params['interp_method'] == None:
        interp_method = ['cubic'] * num_scapes
    else:
        interp_method = params['interp_method']


    dims = params['dims']

    n_rand_pts = params['n_rand_pts']



    return Landscape_stack([random_surface(dims, n_rand_pts, interp_method = interp_method[n], num_hab_types = num_hab_types) for n in range(num_scapes)])




