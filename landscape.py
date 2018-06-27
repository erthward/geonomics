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
import viz

import numpy as np
import numpy.random as r
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter as C
from operator import itemgetter as ig
from scipy import interpolate
from shapely import geometry as g


# ------------------------------------
# CLASSES ---------------------------
# ------------------------------------


class Landscape:
    def __init__(self, dims, raster):
        self.dims = dims
        self.raster = raster
        self.island_mask = False
        self.mask_island_vals = False
        assert type(self.dims) in [tuple, list], "dims must be expressed on a tuple or a list"
        assert type(self.raster) == np.ndarray, "raster should be a numpy.ndarray"

    #####################
    ### OTHER METHODS ###
    #####################

    def show(self, colorbar=True, im_interp_method='nearest', cmap = 'terrain', x=None, y=None, zoom_width=None):
        if self.mask_island_vals:
            mask_val = 1e-7
        else:
            mask_val = None
        plt_lims = viz.get_plt_lims(self, z, y, zoom_width)
        viz.show_rasters(self, colorbar = colorbar, im_interp_method = im_interp_method, cmap = cmap, plt_lims = plt_lims, mask_val = mask_val)


class Landscape_Stack:
    def __init__(self, raster_list, params=None):
        self.scapes = dict(zip(range(len(raster_list)), raster_list))
        assert False not in [raster.__class__.__name__ == 'Landscape' for raster in
                             raster_list], 'All landscapes supplied in raster_list must be of type landscape.Landscape.'
        self.n_scapes = len(raster_list)
        assert len(set([land.dims for land in list(self.scapes.values())])) == 1, 'Dimensions of all landscapes must be even.'
        self.dims = list(self.scapes.values())[0].dims
        self.dim_om = max([len(str(d)) for d in self.dims]) #get the order of magnitude of the larger land dimension (used when zero-padding cell-coorindate strings)
        self.movement_surf = None
        self.movement_surf_approx_len = None
        self.movement_surf_scape_num = None
        self.density_grid_stack = None
        self.n_island_mask_scape = None
        #self.mating_grid = mating_grid.MatingGrid(params=params)


    #####################
    ### OTHER METHODS ###
    #####################

    def show(self, scape_num=None, colorbar=True, cmap='terrain', im_interp_method='nearest', x=None, y=None, zoom_width=None):
        if True in [scape.mask_island_vals for scape in self.scapes.values()]:
            mask_val = 1e-7
        else:
            mask_val = None
        plt_lims = viz.get_plt_lims(land, x, y, zoom_width)
        viz.show_rasters(self, scape_num = scape_num, colorbar = colorbar, im_interp_method = im_interp_method, cmap = cmap, mask_val = mask_val, plt_lims = plt_lims)


    #method for plotting the movement surface (in various formats)
    def show_movement_surf(self, style, x, y, zoom_width=8, scale_fact=4.5, color='black', colorbar = True):
        '''
        'style' can take values: 'hist', 'circ_hist', 'vect', or 'circ_draws'
        '''
        if self.movement_surf is None:
            print('This landscape stack appears to have no movement surface layer. Function not valid.')
            return
        elif style not in ['hist', 'circ_hist', 'vect', 'circ_draws']:
            print("The 'style' argument must take one of the following values: 'hist', 'circ_hist', 'vect', 'circ_draws'")
            return
        elif style == 'hist':
            plt.hist(r.choice(self.movement_surf[i,j,:], size = 10000, replace = True), bins=100, density=True, alpha=0.5)

        else:
            #display the movement-surface raster
            scape_num = self.movement_surf_scape_num
            self.scapes[scape_num].show(zoom_width = zoom_width, x = x, y = y)

            if style == 'circ_hist':
                v, a = np.histogram(r.choice(self.movement_surf[y,x,:], replace = True, size = 7500), bins=15)
                v = v / float(v.sum())
                a = [(a[n] + a[n + 1]) / 2 for n in range(len(a) - 1)]
                xs = [np.cos(a[n]) * 0.5 for n in range(len(a))]
                ys = [np.sin(a[n]) * 0.5 for n in range(len(a))]
                xs = np.array(xs) * v * scale_fact
                ys = np.array(ys) * v * scale_fact
                [plt.plot((x, (x + xs[n])), (y, (y + ys[n])), linewidth=2, color=color) for n in range(len(xs))]
           
            elif style == 'circ_draws':
                pts = [(np.cos(a), np.sin(a)) for a in r.choice(self.movement_surf[y,x,:], size = 1000, replace = True)]
                plt.scatter([pt[0] * 0.5 + x for pt in pts], [pt[1] * 0.5 + y for pt in pts], color='red', alpha=0.1, marker = '.')

            elif style == 'vect':
                def plot_one_cell(i, j):
                    # draw sample of angles from the Gaussian KDE representing the von mises mixture distribution (KDE)
                    samp = self.movement_surf[i,j,:]
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
                [[plot_one_cell(i, j) for i in range(self.movement_surf.shape[0])] for j in range(self.movement_surf.shape[1])]
           

    # method for pickling a landscape stack
    def pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)


class Density_Grid:
    def __init__(self, dims, dim_om, window_width, gi, gj, cells, areas, x_edge, y_edge):

        #same as land.dims
        self.dims = dims

        #order of magnitude of the largest dimension in self.dims (used for zero-padding the cell strings)
        self.dim_om = dim_om

       #window width to be used for grid-cell windows within which population will be counted
        self.window_width = window_width
       
        #self.x_edge/y_edge are True if this grid has cells centered on the edge (i.e. 0 and the dims val) for this dimension, else False
        self.x_edge = x_edge
        self.y_edge = y_edge
       

        #meshgrid arrays of the i and j points associated with each window
        self.gi = gi
        self.gj = gj

        #create array of the i,j grid-cell centerpoints
        self.grid_coords = np.array(list(zip(self.gi.flatten(), self.gj.flatten())))

        #create attriibutes of the grid cells and their areas
        self.cells = cells
        self.areas = areas

        #save itemgetter as an attribute (returns the counted number of instances of each cell from the collections.Counter() dict)
        self.grid_counter = ig(*self.cells) 


    def calc_density(self, x, y):

        #determine the x and y cells within which each individual's x and y coordinates fall
        #(the self.x_edge and self.y_edge corrections determine the correct cells whether the grid includes cells
        #centered on the landscape edges , i.e. self.x_edge or y_edge = True, or not)
        x_cells = (x - self.x_edge*self.window_width/2.)//self.window_width + self.x_edge
        y_cells = (y - self.y_edge*self.window_width/2.)//self.window_width + self.y_edge
        
        #get cell strings from indivudals' cells
        cells = get_cell_strings(y_cells, x_cells, self.dim_om)
        #and get Counter dict of the cell strings
        counts = C(cells)

        #use the itemgetter attribute to get the count for each cell
        grid_counts = self.grid_counter(counts)
        #reshape them into an ndarray 
        grid_counts = np.reshape(grid_counts, self.gi.shape)
        #and divide the array values by the appropriate grid-cell areas to get the densities
        grid_dens = grid_counts/self.areas

        return(grid_dens)



class Density_Grid_Stack:
    def __init__(self, land, window_width = None):
    
        self.dims = land.dims

        #NOTE: This is a quick hack to find the closest whole-number factor of the larger of the landscape
        #dimensions and set to 1/10th of that dimension, and set that as the window_width.
        #This is so the windows fit evenly across the landscape. Otherwise, strange numerical errors are
        #occurring. But when I get the time, I need to puzzle through the density_grid_stack code again and
        #figure out why it throws errors in those cases, and how to rectify it
        if window_width == None:
            facts = [(i,abs(i-0.1*max(self.dims))) for i in range(1,max(self.dims)+1) if max(self.dims)%i == 0]
            closest_fact = min([f[1] for f in facts])
            window_width = [f[0] for f in facts if f[1] == closest_fact][0]

        #set window width of the grid-cell windows within which to count the population
        self.window_width = window_width

        #get the order of magnitude of the larger of the two land dimensions (used to zero-pad the cell-coordinate strings)
        self.dim_om = land.dim_om

        #get meshgrids of the i and j cell-center coordinates of the landscape-raster cells (to be interpolated to for density calculation)
        self.land_gi, self.land_gj = np.meshgrid(np.arange(0, self.dims[0])+0.5, np.arange(0, self.dims[1])+0.5)

        #create inner and outer density grids from the land and window-width
        self.grids = dict([(n,g) for n,g in enumerate(make_density_grids(land, self.window_width))])


    def calc_density(self, x, y):
        #get a concatenated list of the grid-cell center coordinates from both density grids
        pts = np.vstack([self.grids[n].grid_coords for n in range(len(self.grids))])
        #and a concatenated list of the densities calculated for both density grids
        vals = np.hstack([self.grids[n].calc_density(x, y).flatten() for n in range(len(self.grids))])

        #then interpolate from those points and values to the centerpoints of all of the land centerpoints
        dens = interpolate.griddata(pts, vals, (self.land_gj, self.land_gi), method = 'cubic')

        return(dens)


class Movement_Surface:
    def __init__(self, movement_surf, scape_num):
        self.movement_surf = movement_surf
        self.scape_num = scape_num

    def get(x, y, z):
        return(self.movement_surf[y, x, z])


# --------------------------------------
# FUNCTIONS ---------------------------
# --------------------------------------


def random_surface(dims, n_rand_pts, interp_method="cubic", island_val=0, num_hab_types=2, dist='beta', alpha=0.05,
                   beta=0.05):  # requires landscape to be square, such that dim = domain = range
    # NOTE: can use "nearest" interpolation to create random patches of habitat (by integer values); can change num_hab_types to > 2 to create a randomly multi-classed landscape
    # NOTE: I guess this could be used for rectangular landscapes, if the square raster is generated using the larger of the two dimensions, and then the resulting array is subsetted to the landscape's dimensions
    # NOTE: This seems to generate decent, believable random rasters!
    # n_rand_pts/dim ratio:
    # ~0.01-0.05 --> broad, simple gradients
    # ~0.05-0.10 --> slightly more complex, occasionally landscape of low regions divided by a major high region (or vice versa)
    # ~0.10-0.50 --> landscape can be broken up into numerous fragmented patches (though sometimes still relatively homogeneous, with one or two small, extremely different patches
    # ~0.50-1.00 --> more highly fragmented version of above
    max_dim = max(dims)
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
    if dims[0] != dims[1]:
        pass  # NOTE: figure out how to get it to use the dims tuple to subset an approriate size if dims not equal

    return Landscape(dims, I)


def defined_surface(dims, pts, vals, interp_method="cubic",
                    num_hab_types=2):  # pts should be provided as n-by-2 Numpy array, vals as a 1-by-n Numpy array

    # NOTE: There seem to be some serious issues with all of this code, because the resulting landscapes are not quite symmetrical; and making small tweaks (e.g. adding 5 to all input points' coordinates) doesn't just tweak the output landscape but instead completely gets ride of it; I have an intuition that it comes from the code that coerces all raster values to 0 <= val <= 1, becuase that doesn't look it does quite what I originally intended for it to do, but I'm not really sure... anyhow, for now it works for my initial testing purposes

    # NOTE: like the random_surface function, this also requires landscape to be square, but I guess this could be used for rectangular landscapes, if the square raster is generated using the larger of the two dimensions, and then the resulting array is subsetted to the landscape's dimensions

    # NOTE: if discrete habitat patches desired, values should still be fed in as proportions, and will then be multipled by num_hab_types to develop habitat class values
    if interp_method == 'nearest':
        vals = vals * (num_hab_types - 1)

    # add 0.5 to all pts, to center them with respect to the raster display, to make intended symmetry actually symmetrical
    # pts = pts + 0.5

    max_dim = max(dims)
    grid_x, grid_y = np.mgrid[1:max_dim:complex("%ij" % max_dim), 1:max_dim:complex("%ij" % max_dim)]
    I = interpolate.griddata(pts, vals, (grid_x, grid_y), method=interp_method)
    if interp_method == 'nearest':  # i.e., if being used to generate random habitat patches...
        I = I.round().astype(float)
    if interp_method == 'cubic':  # transform to constrain all values to 0 <= val <= 1
        I = I + abs(I.min()) + (
                0.01 * r.rand())  # NOTE: adding a small jitter to keep values from reaching == 0 or == 1,
        # as would likely be the case with linear interpolation
        I = I / (I.max() + (0.01 * r.rand()))
    if dims[0] != dims[1]:
        pass  # NOTE: figure out how to get it to use the dims tuple to subset an approriate size if dims not equal
    return Landscape(dims, I)


def build_scape_stack(params, num_hab_types=2):
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

    dims = params['land']['dims']

    # create rasters for random landscape, if params['land']['rand_land'] == True
    if params['land']['rand_land']:

        n_rand_pts = params['land']['n_rand_pts']

        # if only a single integer provided for n_rand_pts, then use that for all landscape layers (i.e. scale of
        # spatial heterogeneity will be roughly equal for all layers); otherwise, a list or tuple of n_rand_pts could
        #  create different scales of heterogeneity for different landscape layers
        if type(n_rand_pts) == int:
            n_rand_pts = [n_rand_pts] * num_scapes

        land = Landscape_Stack( [random_surface(dims, n_rand_pts[n], interp_method=interp_method[n], num_hab_types=num_hab_types) for n in range(num_scapes)], params=params)

    # or create rasters for defined landscape, if params['land']['rand_land'] == False
    elif not params['land']['rand_land']:

        # get pts
        pts = params['land']['landscape_pt_coords']
        # if only a single array of pts provided, multiply them into a list to use for each separate landscape layer
        if type(pts) == np.ndarray:
            pts = [pts] * num_scapes

        # get vals
        vals = params['land']['landscape_pt_vals']

        land = Landscape_Stack(
            [defined_surface(dims, pts[n], vals[n], interp_method=interp_method[n], num_hab_types=num_hab_types) for n
             in range(num_scapes)])

    # if params['islands']['island_val'] > 0, then use the movement_surf raster to create T/F mask-raster, added as an
    # additional landscape, to be used to kill all individuals straying off 'islands'

    #create the land.density_grid_stack object
    if params['land']['density_grid_window_width'] != None:
        land.density_grid_stack = Density_Grid_Stack(land, window_width = params['land']['density_grid_window_width'])
    else:
        land.density_grid_stack = Density_Grid_Stack(land)

    # create a movement surface, if params call for it
    # NOTE: THIS WILL TAKE A WHILE TO COMPUTER UP-FRONT!
    if params['land']['movement_surf']['movement_surf']:

        land.movement_surf_scape_num = params['land']['movement_surf']['movement_surf_scape_num']
        land.movement_surf_approx_len = params['land']['movement_surf']['movement_surf_approx_len']

        if params['land']['islands']['islands']:

            if params['land']['islands']['island_val'] > 0:
                iv = params['land']['islands']['island_val']

                # zero out the appropriate parts of the movement raster
                movement_rast = land.scapes[land.movement_surf_scape_num].raster
                # zero out the appropriate parts of the movement raster (but not exactly 0, because
                # creates division-by-zero problems in the pop_dynamics calculations)
                movement_rast[movement_rast < iv] = 1e-8


            elif 'island_mask' in params['land']['islands'].keys() and params['land']['islands']['island_mask'] is not None:
                im_file = params['land']['islands']['island_mask']
                ma = np.loads(im_file)
                assert type(
                    ma) == np.ndarray, "The pickled file located at the path provided in params['land']['islands']['island_mask'] does not appear to be a numpy ndarray. "

                # get the movement raster scape number
                movement_rast = land.scapes[land.movement_surf_scape_num].raster
                # zero out the appropriate parts of the movement raster (but not exactly 0, because
                # creates division-by-zero problems in the pop_dynamics calculations)
                movement_rast[ma == 0] = 1e-8

            # replace the movement_surf_scape raster with the updated raster with all outside-habitat values set to 1e-8
            land.scapes[land.movement_surf_scape_num].raster = movement_rast
            # set the Landscape.mask_island_vals flag on the movement surf to True, for plotting purposes
            land.scapes[land.movement_surf_scape_num].mask_island_vals = True

            # then create an island mask and add as the last landscape (to use to quickly cull individuals outside
            # the 'islands')
            island_mask = create_island_mask(dims, land.scapes[land.movement_surf_scape_num].raster)
            island_mask.island_mask = True  # set Lanscape.island_mask attribute to True
            land.scapes[land.n_scapes] = island_mask  # set as the last scape
            land.n_island_mask_scape = land.n_scapes  # the Landscape_Stack n_island_mask_scape attribute to last scape num
            land.n_scapes += 1  # then increment total scape nums by 1

        # create the movement surface, and set it as the land.movement_surf attribute
        import movement
        land.movement_surf = movement.create_movement_surface(land, approx_len = land.movement_surf_approx_len)

    return land


def create_island_mask(dims, scape):
    mask = np.ma.masked_less_equal(scape, 1e-8).mask
    return Landscape(dims, mask)


# function for reading in a pickled landscape stack
def load_pickled_scape_stack(filename):
    import cPickle
    with open(filename, 'rb') as f:
        land = cPickle.load(f)

    return land


####functions for creating density-calculation grids

#create strings from input cell coordinates
def get_cell_strings(gi, gj, dim_om):
    #get strings for both i and j cooridnates, zfilling to the correct order-of-magnitude (of the larger landscape dimension)
    i_strs = [str(int(i)).zfill(dim_om) for i in gi.flatten()]
    j_strs = [str(int(j)).zfill(dim_om) for j in gj.flatten()]

    #join those strings to one single string, unique for each cell
    cells = [''.join(c) for c in list(zip(i_strs,j_strs))]

    return(cells)




#make a density grid, based on the Landscape_Stack, the chosen window-width, 
#and the Boolean arguments dictating whether or not the grid's 
#x- and y-dimension cells should be centered on the landscape 
#edges (i.e. 0 and dims[_])
def make_density_grid(land, ww, x_edge, y_edge):
    
    #half-window width
    hww = ww/2.

    #get land dimensions
    dims = land.dims
    dim_om = land.dim_om
   
    #create a dictionary of cell ranges, one for when cells center on edge values 
    #(i.e. 0 and dims[n] for either dimension), 
    #the other for when they don't 
    #(i.e. run from hww to dims[n] - hww)
    edge_range_dict = {True:  np.arange(0, dims[0]+ww, ww), 
                       False: np.arange(0+hww, dims[0]+hww, ww)}
    
    #create the meshgrid of the centerpoints of neighborhoods (or cells) within which population will be counted
    #(x_edge and y_edge arguments determine whether this grid's 
    #x and y cells are centered on the lanscape edges or not)
    #NOTE: these are expressed as points in continuous space from 0 to each landscape dimension, NOT as cell
    #numbers (which will be calculated below)
    gj, gi = np.meshgrid(edge_range_dict[x_edge], edge_range_dict[y_edge])
    
    #and get flattened lists of the grid's i and j coordinates 
    j = gj.flatten()
    i = gi.flatten()


    #create a single, large Polygon object of the landscape quadrilateral
    land_poly_coords = ((0,0), (dims[0], 0), (dims[0], dims[1]), (0, dims[1]))
    land_poly = g.Polygon(land_poly_coords)
   
    
    #create a list of quadrilaterals centered on each of the points in the grid
    polys = [g.Polygon(((j[n]-hww, i[n]-hww), (j[n]-hww, i[n]+hww), (j[n]+hww, i[n]+hww), (j[n]+hww, i[n]-hww))) for n in range(len(j))]

    #use the Polygons' intersection method to calculate the total area within the landscape that is covered
    #by each grid cell (which will be used as the denominator for calculating neighborhood population densities
    #from neighborhood population counts)
    areas = np.reshape([p.intersection(land_poly).area for p in polys], gj.shape)

    #get lists of the integer-identifiers (i,j in the proper matrix sense) of the cells in each of the grids
    #(so that these can be matched up with individuals' cell numbers when calculating population density while
    #the model is running)
    #(e.g. a cell centered at point y = 0, x = 5 would be cell i = 0,j = 6 if the window-width is 1, 
    #or cell i =0,j = 1 if the window width is 5)
    i_cells = (i - (hww * (y_edge))) // ww + (y_edge)
    j_cells = (j - (hww * (x_edge))) // ww + (x_edge)


    #turn the cell-number integers into cell strings
    #(NOTE: the previous algorithm used to calculate population density was more or less the same, but found
    #individuals' cells by checking numerically whether they were within each window, but was much slower and
    #scaled poorly with decreasing window width, increasing landscape size, and increasing population size.
    #This approach instead uses the floor-divide // on individuals' x and y coordinates to generate
    #cell-number integers for them, converts those to strings, then counts the number of instances of each
    #string for the population and uses those values as the counts of individuals within each density-grid
    #cell, obviating the need to loop across grid dimensions. This performs better and scales much better.
    cells = get_cell_strings(i_cells, j_cells, dim_om)

    #use the above-created data structures to create two Density_Grid objects (which will inhere to the
    #Landscape_Stack as attributes)
    grid = Density_Grid(dims, dim_om, ww, gi, gj, cells, areas, x_edge= x_edge, y_edge = y_edge)

    return(grid)

#create 4 density grids, one for each offset (i.e. each combination of offset by 0 and by 0.5*window_width)
def make_density_grids(land, ww):

    #make a grid for each combo of Booleans for x_edge and y_edge
    g1 = make_density_grid(land, ww, x_edge = True, y_edge = True)
    g2 = make_density_grid(land, ww, x_edge = False, y_edge = False)
    g3 = make_density_grid(land, ww, x_edge = True, y_edge = False)
    g4 = make_density_grid(land, ww, x_edge = False, y_edge = True)

    return(g1, g2, g3, g4)

