#!/usr/bin/python
# landscape.py

'''
##########################################

Module name:          landscape


Module contains:
                      - definition of the Scape and Land classes
                      - function for creating a random landscape, based on input parameters
                      - associated functions


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           12-28-15
Documentation:        URL


##########################################
'''
#geonomics imports
from utils import viz, io, spatial as spt

#other imports
import numpy as np
import numpy.random as r
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter as C
from operator import itemgetter as ig
from scipy import interpolate
from shapely import geometry as g


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

class Scape:
    def __init__(self, rast, name, dim, res=(1,1), ulc = (0,0)):
        self.idx = None
        self.name = name
        self.dim = dim
        assert type(self.dim) in [tuple, list], "dim must be expressed on a tuple or a list"
        #set the resoultion (res; i.e. cell-size) and upper-left corner (ulc) to defaults; will be reset if landscape read in from a GIS raster
        self.res = res  
        self.ulc = ulc
        self.rast = rast
        assert type(self.rast) == np.ndarray, "rast should be a numpy.ndarray"
        self.scale_min = None
        self.scale_max = None
        self.island_mask = False
        self.mask_island_vals = False

    #####################
    ### OTHER METHODS ###
    #####################

    def show(self, colorbar=True, im_interp_method='nearest', cmap = 'terrain', x=None, y=None, zoom_width=None, vmin=None, vmax=None):
        if self.mask_island_vals:
            mask_val = 1e-7
        else:
            mask_val = None
        plt_lims = viz.get_plt_lims(self, x, y, zoom_width)
        viz.show_rasters(self, colorbar = colorbar, im_interp_method = im_interp_method, cmap = cmap, plt_lims = plt_lims, mask_val = mask_val, vmin = vmin, vmax = vmax)


class Land(dict):

    #######################
    ### SPECIAL METHODS ###
    #######################

    def __init__(self, scapes, params=None, res=(1,1), ulc = (0,0)):
        #check the scapes dict is correct, then update the Land with it
        assert False not in [scape.__class__.__name__ == 'Scape' for scape in scapes.values()], 'ERROR: All landscapes supplied in scapes must be of type landscape.Scape.'
        self.update(scapes)
        #set the scapes attribute to the dict values (it still changes dynamically as things update, it runs
        #3x faster than calling values() during the model, and it makes code nicer looking
        self.scapes = self.values()
        #set the number of scapes in the stack
        self.n_scapes = len(self)
        #set the idx attributes on the scapes
        [setattr(v, 'idx', k) for k, v in self.items()]

        #check that scapes' dimensions are all the same
        assert len(set([land.dim for land in list(self.values())])) == 1, 'ERROR: Dimensions of all landscapes must be equal.'
        #then set the dim to the 0th landscape's dim 
        self.dim = list(self.values())[0].dim
        #get the order of magnitude of the larger land dimension (used when zero-padding cell-coorindate strings)
        self.dim_om = max([len(str(d)) for d in self.dim]) 
        #set the resoultion (res; i.e. cell-size) and upper-left corner (ulc) to 0th scape's values (equality
        #across scapes should already have been checked); will be reset if landscape read in from a GIS raster
        self.res = [*self.values()][0].res
        self.ulc = [*self.values()][0].ulc

        #create a changer attribute that defaults to None but will be set to an ops.change.Land_Changer object
        #if params call for it
        self.changer = None

        self.island_mask_scape_num = None

    #define the __str__ and __repr__ special methods
    #NOTE: this doesn't excellently fit the Python docs' specification for __repr__; I should massage this #some more once I'm done writing the codebase
    def __str__(self):
        type_str = str(type(self))
        scapes_str = '\n%i Scapes:\n' % self.n_scapes
        scapes_str = scapes_str + ',\n'.join(sorted(['\t' + str(k) + ': ' + str(v) for k,v in self.items()])) + '},'
        vars_str = "\nParameters:\n\t" + ',\n\t'.join(sorted([str(k) + ': ' + str(v) for k,v in vars(self).items()])) + '}}'
        return '\n'.join([type_str, scapes_str, vars_str])

    def __repr__(self):
        repr_str = self.__str__()
        return repr_str

    #####################
    ### OTHER METHODS ###
    #####################

    #method to set a raster
    def set_raster(self, scape_num, rast):
        self[scape_num].rast = rast

    #method to make landscape changes
    def make_change(self, t):
        self.changer.make_change(t)

    #method to plot the landscape (or just a certain scape)
    def show(self, scape_num=None, colorbar=True, cmap='terrain', im_interp_method='nearest', x=None, y=None, zoom_width=None, vmin=None, vmax=None):
        if True in [scape.mask_island_vals for scape in self.values()]:
            mask_val = 1e-7
        else:
            mask_val = None
        plt_lims = viz.get_plt_lims(self, x, y, zoom_width)
        viz.show_rasters(self, scape_num = scape_num, colorbar = colorbar, im_interp_method = im_interp_method, cmap = cmap, mask_val = mask_val, plt_lims = plt_lims, vmin = vmin, vmax = vmax)


    #method for plotting the movement surface (in various formats)
    def show_movement_surface(self, style, x, y, zoom_width=8, scale_fact=4.5, color='black', colorbar = True):
        '''
        'style' can take values: 'hist', 'circ_hist', 'vect', or 'circ_draws'
        '''
        if self.move_surf is None:
            print('ERROR: This landscape stack appears to have no movement surface layer. Function not valid.')
            return
        elif style not in ['hist', 'circ_hist', 'vect', 'circ_draws']:
            print("ERROR: The 'style' argument must take one of the following values: 'hist', 'circ_hist', 'vect', 'circ_draws'")
            return
        elif style == 'hist':
            plt.hist(r.choice(self.move_surf[i,j,:], size = 10000, replace = True), bins=100, density=True, alpha=0.5)

        else:
            #display the movement-surface raster
            scape_num = self.move_surf_scape_num
            self[scape_num].show(zoom_width = zoom_width, x = x, y = y)

            if style == 'circ_hist':
                v, a = np.histogram(r.choice(self.move_surf[y,x,:], replace = True, size = 7500), bins=15)
                v = v / float(v.sum())
                a = [(a[n] + a[n + 1]) / 2 for n in range(len(a) - 1)]
                xs = [np.cos(a[n]) * 0.5 for n in range(len(a))]
                ys = [np.sin(a[n]) * 0.5 for n in range(len(a))]
                xs = np.array(xs) * v * scale_fact
                ys = np.array(ys) * v * scale_fact
                [plt.plot((x, (x + xs[n])), (y, (y + ys[n])), linewidth=2, color=color) for n in range(len(xs))]
           
            elif style == 'circ_draws':
                pts = [(np.cos(a), np.sin(a)) for a in r.choice(self.move_surf[y,x,:], size = 1000, replace = True)]
                plt.scatter([pt[0] * 0.5 + x for pt in pts], [pt[1] * 0.5 + y for pt in pts], color='red', alpha=0.1, marker = '.')

            elif style == 'vect':
                def plot_one_cell(i, j):
                    # draw sample of angles from the Gaussian KDE representing the von mises mixture distribution (KDE)
                    samp = self.move_surf[i,j,:]
                    # create lists of the x and y (i.e. cos and sin) components of each angle in the sample
                    x_vects = np.cos(samp)
                    y_vects = np.sin(samp)
                    # define the dx and dy distances used to the position the arrowhead
                    # (divide by sqrt(2)/2, to scale to the size of half of the diagonal of a cell)
                    dx = np.mean(x_vects) / np.sqrt(2)
                    dy = np.mean(y_vects) / np.sqrt(2)
                    # now plot the arrow
                    plt.arrow(j, i, dx, dy, alpha=0.75, color='black', head_width=0.24, head_length=0.32)

                # call the internally defined function as a nested list comprehension for all raster cells, which I believe should do its best to vectorize the whole operation
                [[plot_one_cell(i, j) for i in range(self.move_surf.shape[0])] for j in range(self.move_surf.shape[1])]
           

    # method for pickling a landscape stack
    def write_pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

def make_random_scape(dim, n_rand_pts, interp_method="cubic", island_val=0, num_hab_types=2, dist='beta', alpha=0.05,
                   beta=0.05):  # requires landscape to be square, such that dim = domain = range
    # NOTE: can use "nearest" interpolation to create random patches of habitat (by integer values); can change num_hab_types to > 2 to create a randomly multi-classed landscape
    # NOTE: I guess this could be used for rectangular landscapes, if the square raster is generated using the larger of the two dimensions, and then the resulting array is subsetted to the landscape's dimensions
    # NOTE: This seems to generate decent, believable random rasters!
    # n_rand_pts/dim ratio:
    # ~0.01-0.05 --> broad, simple gradients
    # ~0.05-0.10 --> slightly more complex, occasionally landscape of low regions divided by a major high region (or vice versa)
    # ~0.10-0.50 --> landscape can be broken up into numerous fragmented patches (though sometimes still relatively homogeneous, with one or two small, extremely different patches
    # ~0.50-1.00 --> more highly fragmented version of above
    max_dim = max(dim)
    if interp_method == 'nearest':
        if dist == 'unif':
            vals = r.rand(n_rand_pts) * (num_hab_types - 1)
        elif dist == 'beta':
            vals = r.beta(alpha, beta, n_rand_pts) * (num_hab_types - 1)
    else:
        if dist == 'unif':
            vals = r.rand(n_rand_pts)
        elif dist == 'beta':
            vals = r.beta(alpha, beta, n_rand_pts)
    pts = r.normal(max_dim / 2, max_dim * 2, [n_rand_pts,
                                              2])  # selects seed points from well outside the eventaul landscape, to ensure interpolation across area of interest
    grid_x, grid_y = np.mgrid[1:max_dim:complex("%ij" % max_dim), 1:max_dim:complex(
        "%ij" % max_dim)]  # by this function's definition, the use of complex numbers in here specifies the number of steps desired
    I = interpolate.griddata(pts, vals, (grid_x, grid_y), method=interp_method)
    if interp_method == 'nearest':  # i.e., if being used to generate random habitat patches...
        I = I.round().astype(float)
    if interp_method == 'cubic':  # transform to constrain all values to 0 <= val <= 1
        I = I + abs(I.min()) + (
                0.01 * r.rand())  # NOTE: adding a small jitter to keep values from reaching == 0 or == 1, as would likely be the case with linear interpolation
        I = I / (I.max() + (0.01 * r.rand()))
    if dim[0] != dim[1]:
        pass  # NOTE: figure out how to get it to use the dim tuple to subset an approriate size if dim not equal

    return I


def make_defined_scape(dim, pts, vals, interp_method="cubic", num_hab_types=2):  # pts should be provided as n-by-2 Numpy array, vals as a 1-by-n Numpy array

    # NOTE: There seem to be some serious issues with all of this code, because the resulting landscapes are not quite symmetrical; and making small tweaks (e.g. adding 5 to all input points' coordinates) doesn't just tweak the output landscape but instead completely gets ride of it; I have an intuition that it comes from the code that coerces all raster values to 0 <= val <= 1, becuase that doesn't look it does quite what I originally intended for it to do, but I'm not really sure... anyhow, for now it works for my initial testing purposes

    # NOTE: like the make_random_scape function, this also requires landscape to be square, but I guess this could be used for rectangular landscapes, if the square raster is generated using the larger of the two dimensions, and then the resulting array is subsetted to the landscape's dimensions

    # NOTE: if discrete habitat patches desired, values should still be fed in as proportions, and will then be multipled by num_hab_types to develop habitat class values
    if interp_method == 'nearest':
        vals = vals * (num_hab_types - 1)

    # add 0.5 to all pts, to center them with respect to the raster display, to make intended symmetry actually symmetrical
    # pts = pts + 0.5

    max_dim = max(dim)
    grid_x, grid_y = np.mgrid[1:max_dim:complex("%ij" % max_dim), 1:max_dim:complex("%ij" % max_dim)]
    I = interpolate.griddata(pts, vals, (grid_x, grid_y), method=interp_method)
    if interp_method == 'nearest':  # i.e., if being used to generate random habitat patches...
        I = I.round().astype(float)
    if interp_method == 'cubic':  # transform to constrain all values to 0 <= val <= 1
        I = I + abs(I.min()) + (
                0.01 * r.rand())  # NOTE: adding a small jitter to keep values from reaching == 0 or == 1,
        # as would likely be the case with linear interpolation
        I = I / (I.max() + (0.01 * r.rand()))
    if dim[0] != dim[1]:
        pass  # NOTE: figure out how to get it to use the dim tuple to subset an approriate size if dim not equal
    return I


def make_land(params, num_hab_types=2):
    # NOTE: If a multi-class (rather than binary) block-habitat raster would be of interest, would need to make
    # num_hab_types customizable)

    # grab necessary parameters from the params dict
    n_scapes = len(params.land.scapes)

    dim = params.land.main.dim
    res = params.land.main.res
    ulc = params.land.main.ulc

    #create a dictionary to hold all the scapes to be created
    scapes = {}

    #create a list to hold the information for any GIS rasters to be created (because they will all be created
    #after the loop over the scapes, so that it's simpler to check agreement among raster resoultions and registrations
    gis_scape_params = {'name': [], 'gis_scape_num':[], 'gis_filepath':[], 'gis_scale_max_val':[], 'gis_scale_min_val':[]}

    #then loop over the scapes in params.land.scapes and create each one
    for n, (k, scape_params) in enumerate(params.land.scapes.items()):

        #get the init parameters
        init_params = scape_params.init

        #get the scape name
        name = init_params.name

        #get the interpolation method (which would be used for a random scape or a scape defined by a points array)
        interp_method = init_paramsinterp_method
        #and set to 'cubic' if None is params
        if interp_method == None:
            interp_method = 'cubic'

        #if this scape is to be randomly interpolated
        if init_params.rand:
            #get the number of random points to interpolate from
            n_rand_pts = init_params.n_rand_pts
            #create the random scape and add it to the scapes dict
            scapes[n] = Scape(make_random_scape(dim, n_rand_pts, interp_method=interp_method, num_hab_types=num_hab_types), name = name, dim = dim, res = res, ulc = ulc)

        # or create rasters for defined landscape, if params['land']['rand_land'] == False
        elif not init_params.rand:
            
            #create a GIS scape, if a file is provided
            if init_params.gis_filepath is not None and init_params.scape_pt_coords is None and init_params.scape_pt_vals is None):
                #get the GIS params
                filepath = init_params.gis_filepath
                scale_max = init_params.gis_scale_max
                scale_min = init_params.gis_scale_min
                #then store the number, key, and values of the layer to be replaced with a GIS raster after the Land is
                #created (NOTE: easier to do this all at once afterward, to check that the resolution and
                #registration of the rasters all agree)
                scapes[n] = None
                gis_scape_params['gis_scape_num'].append(n)
                gis_scape_params['name'].append(name)
                [gis_scape_params[k].append(init_params[k]) for k in init_params.keys() if k.startswith('gis_')]

            #or create a scape defined by n pts and n vals, if provided 
            elif init_params.defined_scape_pt_coords is not None and init_params.defined_scape_pt_vals is not None and init_params.gis_filepath is None:
                # get pts
                pts = init_params.defined_scape_pt_coords
                # get vals
                vals = init_params.defined_scape_pt_vals
                #create the defined raster
                scapes[n] = Scape(make_defined_scape(dim, pts, vals, interp_method=interp_method, num_hab_types=num_hab_types), name = name, dim = dim, res = res, ulc = ulc)

            #or if some parameters are provided for both a defined scape and a GIS scape, return an informative error
            elif (init_params.defined_scape_pt_coords is not None or init_params.defined_scape_pt_vals is not None) and init_params.gis_filepath is not None:
                raise ValueError('ERROR: For the %ith scape (key: %s) some parameters where provided for both a GIS scape ("gis_<>" parameters) and a point-defined scape ("defined_scape_<>" parameters).  Please provide parameters for only one or the other.' % (n, str(k)))

            #or if no parameters are provided for a defined scape or a GIS scape, return an informative error 
            #(and rand was already checked False, above) return an informative error
            else:
                raise ValueError("ERROR: From the parameters provided in params['land']['scapes'] it could not be determined whether the %ith scape (key: %s) should be randomly interpolated (a 'random scape'), interpolated from a defined set of valued points (a 'point-defined scape'), or read in from a GIS raster layer (a 'gis scape'). Configure the parameters for this scape to one of these options." % (n, str(k))) 

    #create the land object
    land = Land(scapes, params=params)

    #then set the necessary layers to their GIS rasters
    if len(gis_scape_params) > 0:
        set_gis_rasters(land, **gis_scape_params)

    #grab the change parameters into a dictionary of scape_num:events:events_params hierarchical
    change_params = {k:v.change for k,v in params.land.scapes.items() if 'change' in v.keys()} 
    #and then replace the land.changer attribute with an ops.change.Land_Changer object, if params call for it
    if len(change_params) > 0:
        land.changer = ops.change.Land_Changer(land, change_params) 

    return land


def set_gis_rasters(land, name, gis_scape_num, gis_filepath, gis_scale_min_val, gis_scale_max_val):
    assert len(gis_scape_num) == len(gis_filepath), 'ERROR: Parameters provide a different number of GIS raster files to read in than of scape numbers to set them to .'
    res_list = []
    ulc_list = []
    rast_arrays = []
    for n,filepath in enumerate(gis_filepath):
        #get the array, dim, top-left, and res from io.read_raster
        rast_array, rast_dim, rast_ulc, rast_res = io.read_raster(filepath)
        #check that the dimensions are right
        assert rast_dim == land.dim, 'ERROR: land.dim and the dimensions of the input raster %s appear to differ. Please clip %s to the correct dimensions and try again. %s has dimensions (%i, %i), but land has dimensions (%i,%i)' % (filepath, filepath, filepath, rast_dim[0], rast_dim[1], land.dim[0], land.dim[1])
        #scale the raster to 0<=x<=1
        rast_array, scale_min, scale_max = spt.scale_raster(rast_array, min_inval = gis_scale_min_val[n], max_inval = gis_scale_max_val[n])
        #add the rast_array to the rast_arrays list
        rast_arrays.append(rast_array)
        #and add the land's ulc and res attributes to the appropriate lists
        res_list.append(tuple(rast_res))
        ulc_list.append(tuple(rast_ulc))
        #and update scale_min_vals and scale_max_vals if the values used in spt.scale_raster() don't equal the
        #values from params (should only occur if the params values provided were None)
        if scale_min != gis_scale_min_val[n]:
            gis_scale_min_val[n] = scale_min
        if scale_max != gis_scale_max_val[n]:
            gis_scale_max_val[n] = scale_max

    #check that there is only one unique tuple each in res_list and ulc list (otherwise rasters are not in the same projection and/or registration)
    res_C = C(res_list)
    ulc_C = C(ulc_list)
    assert len(res_C) == 1 and list(res_C.values())[0] == len(res_list), 'ERROR: Some of the GIS rasters read in appear to differ in resolution.'
    assert len(ulc_C) == 1 and list(ulc_C.values())[0] == len(ulc_list), 'ERROR: Some of the GIS rasters read in appear to differ in registration (i.e. top-left corner).'
    #assign the land.res and land.ulc attributes, and the scapes' identical attributes, accordingly
    land.res = res_list[0]
    land.ulc = ulc_list[0]
    #create and set the scapes 
    for scape_num, name, rast in zip(gis_scape_num, name, rast_arrays):
        land[scape_num] = Scape(rast, name = name, dim = land.dim, res = land.res, ulc = land.ulc)
    [setattr(scape, 'res', land.res) for scape in land.values()]
    [setattr(scape, 'ulc', land.ulc) for scape in land.values()]
    [setattr(land[n], 'scale_min', scale_min) for n, scale_min in zip(gis_scape_num, scale_min_vals)]
    [setattr(land[n], 'scale_max', scale_max) for n, scale_max in zip(gis_scape_num, scale_max_vals)]

    
# function for reading in a pickled landscape stack
def read_pickled_land(filename):
    import cPickle
    with open(filename, 'rb') as f:
        land = cPickle.load(f)
    return land

