#!/usr/bin/python
# landscape.py

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
#geonomics imports
from utils import viz, io
from utils import spatial as spt

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

class Landscape:
    def __init__(self, rast, dim, res=(1,1), ulc = (0,0)):
        self.idx = None
        self.dim = dim
        #set the resoultion (res; i.e. cell-size) and upper-left corner (ulc) to defaults; will be reset if landscape read in from a GIS raster
        self.res = res  
        self.ulc = ulc
        self.rast = rast
        self.scale_min = None
        self.scale_max = None
        self.island_mask = False
        self.mask_island_vals = False
        assert type(self.dim) in [tuple, list], "dim must be expressed on a tuple or a list"
        assert type(self.rast) == np.ndarray, "rast should be a numpy.ndarray"

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


class Landscape_Stack:
    def __init__(self, array_list, params=None, res=(1,1), ulc = (0,0)):
        self.scapes = dict(zip(range(len(array_list)), array_list))
        assert False not in [scape.__class__.__name__ == 'Landscape' for scape in
                             array_list], 'All landscapes supplied in array_list must be of type landscape.Landscape.'
        #set the scapes' idx attributes
        [setattr(v, 'idx', k) for k, v in self.scapes.items()]
        self.n_scapes = len(array_list)
        assert len(set([land.dim for land in list(self.scapes.values())])) == 1, 'Dimensions of all landscapes must be even.'
        #set the dim to the 0th landscape's dim (should already have been checked that all scapes' dim attributes are equal)
        self.dim = list(self.scapes.values())[0].dim
        #get the order of magnitude of the larger land dimension (used when zero-padding cell-coorindate strings)
        self.dim_om = max([len(str(d)) for d in self.dim]) 
        #set the resoultion (res; i.e. cell-size) and upper-left corner (ulc) to 0th scape's values (equality
        #across scapes should already have been checked); will be reset if landscape read in from a GIS raster
        self.res = list(self.scapes.values())[0].res
        self.ulc = list(self.scapes.values())[0].ulc
        self.move_surf = None
        self.move_surf_approx_len = None
        self.move_surf_scape_num = None
        self.density_grid_stack = None
        self.n_island_mask_scape = None
        #self.mating_grid = mating_grid.MatingGrid(params=params)


    #####################
    ### OTHER METHODS ###
    #####################

    def set_raster(self, scape_num, rast):
        self.scapes[scape_num].rast = rast

    def show(self, scape_num=None, colorbar=True, cmap='terrain', im_interp_method='nearest', x=None, y=None, zoom_width=None, vmin=None, vmax=None):
        if True in [scape.mask_island_vals for scape in self.scapes.values()]:
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
            print('This landscape stack appears to have no movement surface layer. Function not valid.')
            return
        elif style not in ['hist', 'circ_hist', 'vect', 'circ_draws']:
            print("The 'style' argument must take one of the following values: 'hist', 'circ_hist', 'vect', 'circ_draws'")
            return
        elif style == 'hist':
            plt.hist(r.choice(self.move_surf[i,j,:], size = 10000, replace = True), bins=100, density=True, alpha=0.5)

        else:
            #display the movement-surface raster
            scape_num = self.move_surf_scape_num
            self.scapes[scape_num].show(zoom_width = zoom_width, x = x, y = y)

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


def make_scape_stack(params, num_hab_types=2):
    # NOTE: If a multi-class (rather than binary) block-habitat raster would be of interest, would need to make
    # num_hab_types customizable)

    # grab necessary parameters from the params dict

    if params['land']['num_scapes'] == None:
        num_scapes = 1
    else:
        num_scapes = params['land']['num_scapes']

    if params['land']['interp_method'] == None:
        interp_method = ['cubic'] * num_scapes
    else:
        interp_method = params['land']['interp_method']

    dim = params['land']['dim']
    res = params['land']['res']
    ulc = params['land']['ulc']

    # create rasters for random landscape, if params['land']['rand_land'] == True
    if params['land']['rand_land']:

        n_rand_pts = params['land']['n_rand_pts']

        # if only a single integer provided for n_rand_pts, then use that for all landscape layers (i.e. scale of
        # spatial heterogeneity will be roughly equal for all layers); otherwise, a list or tuple of n_rand_pts could
        #  create different scales of heterogeneity for different landscape layers
        if type(n_rand_pts) == int:
            n_rand_pts = [n_rand_pts] * num_scapes
        
        scape_list = [Landscape(make_random_scape(dim, n_rand_pts[n], interp_method=interp_method[n], num_hab_types=num_hab_types), dim = dim, res = res, ulc = ulc) for n in range(num_scapes)]
        land = Landscape_Stack(scape_list, params=params)

    # or create rasters for defined landscape, if params['land']['rand_land'] == False
    elif not params['land']['rand_land']:

        # get pts
        pts = params['land']['landscape_pt_coords']
        # if only a single array of pts provided, multiply them into a list to use for each separate landscape layer
        if type(pts) == np.ndarray:
            pts = [pts] * num_scapes

        # get vals
        vals = params['land']['landscape_pt_vals']

        scape_list = [Landscape(make_defined_scape(dim, pts[n], vals[n], interp_method=interp_method[n], num_hab_types=num_hab_types), dim = dim, res = res, ulc = ulc) for n in range(num_scapes)]
        land = Landscape_Stack(scape_list, params = params)

    #replace certain scapes with GIS layers, and reset the scapes' ulc and res attributes, if necessary
    gis_params = params['land']['gis']
    if len(gis_params['scape_nums']) == len(gis_params['filepaths']) and len(gis_params['scape_nums']) > 0 and list(C(gis_params['scape_nums']).keys()) != [None] and list(C(gis_params['filepaths']).keys()) != [None]:
        set_gis_rasters(land, **gis_params)

    # if params['islands']['island_val'] > 0, then use the move_surf raster to create T/F mask-raster, added as an
    # additional landscape, to be used to kill all individuals straying off 'islands'

    #create the land.density_grid_stack object
    if params['land']['density_grid_window_width'] != None:
        land.density_grid_stack = spt.Density_Grid_Stack(land, window_width = params['land']['density_grid_window_width'])
    else:
        land.density_grid_stack = spt.Density_Grid_Stack(land)

    # create a movement surface, if params call for it
    # NOTE: THIS WILL TAKE A WHILE TO COMPUTER UP-FRONT!
    if params['land']['move_surf']['move_surf']:

        land.move_surf_scape_num = params['land']['move_surf']['move_surf_scape_num']
        land.move_surf_approx_len = params['land']['move_surf']['move_surf_approx_len']

        if params['land']['islands']['islands']:

            if params['land']['islands']['island_val'] > 0:
                iv = params['land']['islands']['island_val']

                # zero out the appropriate parts of the movement raster
                movement_rast = land.scapes[land.move_surf_scape_num].rast
                # zero out the appropriate parts of the movement raster (but not exactly 0, because
                # creates division-by-zero problems in the pop_dynamics calculations)
                movement_rast[movement_rast < iv] = 1e-8


            elif 'island_mask' in params['land']['islands'].keys() and params['land']['islands']['island_mask'] is not None:
                im_file = params['land']['islands']['island_mask']
                ma = np.loads(im_file)
                assert type(
                    ma) == np.ndarray, "The pickled file located at the path provided in params['land']['islands']['island_mask'] does not appear to be a numpy ndarray. "

                # get the movement raster scape number
                movement_rast = land.scapes[land.move_surf_scape_num].rast
                # zero out the appropriate parts of the movement raster (but not exactly 0, because
                # creates division-by-zero problems in the pop_dynamics calculations)
                movement_rast[ma == 0] = 1e-8

            # replace the move_surf_scape raster with the updated raster with all outside-habitat values set to 1e-8
            land.scapes[land.move_surf_scape_num].rast = movement_rast
            # set the Landscape.mask_island_vals flag on the movement surf to True, for plotting purposes
            land.scapes[land.move_surf_scape_num].mask_island_vals = True

            # then create an island mask and add as the last landscape (to use to quickly cull individuals outside
            # the 'islands')
            island_mask = create_island_mask(dim, land.scapes[land.move_surf_scape_num].rast)
            island_mask.island_mask = True  # set Lanscape.island_mask attribute to True
            land.scapes[land.n_scapes] = island_mask  # set as the last scape
            land.n_island_mask_scape = land.n_scapes  # the Landscape_Stack n_island_mask_scape attribute to last scape num
            land.n_scapes += 1  # then increment total scape nums by 1

        # create the movement surface, and set it as the land.move_surf attribute
        from ops import movement
        land.move_surf = spt.Movement_Surface(land, kappa = 12)

    return land


def make_island_mask(scape):
    mask = np.ma.masked_less_equal(scape, 1e-8).mask
    return Landscape(mask)


def set_gis_rasters(land, scape_nums, filepaths, scale_min_vals, scale_max_vals):
    #make scape_num and filepath lists if they're not
    if type(scape_nums) is int:
        scape_nums = [scape_nums]
    if type(filepaths) is str:
        filepaths = [filepaths]
    assert len(scape_nums) == len(filepaths), 'ERROR: Different numbers of GIS raster file numbers and scape numbers to replace with them have been provided.'
    res_list = []
    ulc_list = []
    rast_arrays = []
    for n,filepath in enumerate(filepaths):
        #get the array, dim, top-left, and res from io.read_raster
        rast_array, rast_dim, rast_ulc, rast_res = io.read_raster(filepath)
        #check that the dimensions are right
        assert rast_dim == land.dim, 'ERROR: land.dim and the dimensions of the input raster %s appear to differ. Please clip %s to the correct dimensions and try again. %s has dimensions (%i, %i), but land has dimensions (%i,%i)' % (filepath, filepath, filepath, rast_dim[0], rast_dim[1], land.dim[0], land.dim[1])
        #scale the raster to 0<=x<=1
        rast_array, scale_min, scale_max = spt.scale_raster(rast_array, min_inval = scale_min_vals[n], max_inval = scale_max_vals[n])
        #add the rast_array to the rast_arrays list
        rast_arrays.append(rast_array)
        #and add the land's ulc and res attributes to the appropriate lists
        res_list.append(tuple(rast_res))
        ulc_list.append(tuple(rast_ulc))
        #and update scale_min_vals and scale_max_vals if the values used in spt.scale_raster() don't equal the
        #values from params (should only occur if the params values provided were None)
        if scale_min != scale_min_vals[n]:
            scale_min_vals[n] = scale_min
        if scale_max != scale_max_vals[n]:
            scale_max_vals[n] = scale_max


    #check that there is only one unique tuple each in res_list and ulc list (otherwise rasters are not in the same projection and/or registration)
    res_C = C(res_list)
    ulc_C = C(ulc_list)
    assert len(res_C) == 1 and list(res_C.values())[0] == len(res_list), 'ERROR: Some of the GIS rasters read in appear to differ in resolution.'
    assert len(ulc_C) == 1 and list(ulc_C.values())[0] == len(ulc_list), 'ERROR: Some of the GIS rasters read in appear to differ in registration (i.e. top-left corner).'
   
    #set the land's correct scape to this raster
    [land.set_raster(scape_num, rast) for scape_num, rast in zip(scape_nums, rast_arrays)]
    #then assign the land.res and land.ulc attributes, and the scapes' identical attributes, accordingly
    land.res = res_list[0]
    land.ulc = ulc_list[0]
    [setattr(scape, 'res', land.res) for scape in land.scapes.values()]
    [setattr(scape, 'ulc', land.ulc) for scape in land.scapes.values()]
    [setattr(land.scapes[n], 'scale_min', scale_min) for n, scale_min in zip(scape_nums, scale_min_vals)]
    [setattr(land.scapes[n], 'scale_max', scale_max) for n, scale_max in zip(scape_nums, scale_max_vals)]

    
# function for reading in a pickled landscape stack
def read_pickled_scape_stack(filename):
    import cPickle
    with open(filename, 'rb') as f:
        land = cPickle.load(f)
    return(land)

