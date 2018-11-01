#!/usr/bin/python
# spatial.py

'''
##########################################

Module name:          utils.spatial


Module contains:
                      - general-use spatial classes and functions, for use
                        in a variety of operations


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           07-04-18
Documentation:        URL


##########################################
'''

#geonomics imports

#other imports
import numpy as np
import numpy.random as r
import random
from numpy import sin, cos, pi
from numpy.random import vonmises as r_vonmises
from scipy.stats import vonmises as s_vonmises
s_vonmises.a = -np.inf
s_vonmises.b = np.inf
from collections import Counter as C
from copy import deepcopy
from operator import itemgetter as ig
from itertools import repeat, starmap
from scipy import interpolate
from scipy.spatial import cKDTree
from shapely import geometry as g
import warnings

#optional imports
try:
    from nlmpy import nlmpy
except Exception as e:
    warnings.warn(("Unable to import module 'nlmpy.nlmpy'. Will not be able to"
        " use it package to make NLMpy rasters (neutral landscape models. The "
        "following error was raised:\n\t%s\n\n") % e)


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

class _DensityGrid:
    def __init__(self, dim, dim_om, window_width, gi, gj, cells, areas,
                                                            x_edge, y_edge):

        #same as land.dim
        self.dim = dim

        #resolution (i.e. cell-size); defaults to 1
        self.res = 1

        #order of magnitude of the largest dimension in self.dim
        #(used for zero-padding the cell strings)
        self.dim_om = dim_om

        #window width to be used for grid-cell windows within which
        #population will be counted
        self.window_width = window_width

        #self.x_edge/y_edge are True if this grid has cells centered
        #on the edge (i.e. 0 and the dim val) for this dimension, else False
        self.x_edge = x_edge
        self.y_edge = y_edge


        #meshgrid arrays of the i and j points associated with each window
        self.gi = gi
        self.gj = gj

        #create array of the i,j grid-cell centerpoints
        self.grid_coords = np.array(list(zip(self.gi.flatten(),
                                                        self.gj.flatten())))

        #create attriibutes of the grid cells and their areas
        self.cells = cells
        self.areas = areas

        #save itemgetter as an attribute (returns the counted
        #number of instances of each cell from the collections.Counter() dict)
        self.grid_counter = ig(*self.cells)


    def _calc_density(self, x, y):
        #determine the x and y cells within which each individual's
        #x and y coordinates fall (the self.x_edge and self.y_edge
        #corrections determine the correct cells whether the grid 
        #includes cells centered on the landscape edges,
        #i.e. self.x_edge or y_edge = True, or not)
        x_cells = (x - self.x_edge*self.window_width/2.
                                            )//self.window_width + self.x_edge
        y_cells = (y - self.y_edge*self.window_width/2.
                                            )//self.window_width + self.y_edge

        #get cell strings from indivudals' cells
        cells = _make_cell_strings(y_cells, x_cells, self.dim_om)
        #and get Counter dict of the cell strings
        counts = C(cells)

        #use the itemgetter attribute to get the count for each cell
        grid_counts = self.grid_counter(counts)
        #reshape them into an ndarray 
        grid_counts = np.reshape(grid_counts, self.gi.shape)
        #and divide the array values by the appropriate
        #grid-cell areas to get the densities
        grid_dens = grid_counts/self.areas

        return grid_dens


class _DensityGridStack:
    def __init__(self, land, window_width = None):

        #dimensions
        self.dim = land.dim

        #resolution (i.e. cell-size)
        self.res = land.res


        #NOTE: This is a quick hack to find the closest whole-number
        #factor of the larger of the land dimensions and set to 1/10th
        #of that dimension, and set that as the window_width.
        #This is so the windows fit evenly across the land. Otherwise,
        #strange numerical errors are occurring. But when I get the 
        #time, I need to puzzle through the density_grid_stack code again and
        #figure out why it throws errors in those cases, and how to rectify it
        if window_width == None:
            facts = [(i,abs(i-0.1*max(self.dim))) for i in range(1,
                                    max(self.dim)+1) if max(self.dim)%i == 0]
            closest_fact = min([f[1] for f in facts])
            window_width = [f[0] for f in facts if f[1] == closest_fact][0]

        #TODO: Should figure out how to make attributes like this one
        #(window_width) unchangeable by the user, because the ability to change
        #them gives the impression that density will be calculated with the
        #newly changed value whereas in fact it will still always be calculated
        #with whatever window_width was initially used to create the grid stack

        #set window width of the grid-cell windows within which
        #to count the population
        self.window_width = window_width

        #get the order of magnitude of the larger of the two land
        #dimensions (used to zero-pad the cell-coordinate strings)
        self.dim_om = land._dim_om

        #get meshgrids of the i and j cell-center coordinates of the
        #landscape-raster cells (to be interpolated to for density calculation)
        self.land_gi, self.land_gj = np.meshgrid(np.arange(0,
                            self.dim[0])+0.5, np.arange(0, self.dim[1])+0.5)

        #create inner and outer density grids from the land and window-width
        self.grids = dict([(n,g) for n,g in enumerate(_make_density_grids(land,
                                                        self.window_width))])


    def _calc_density(self, x, y):
        #get a concatenated list of the grid-cell center coordinates
        #from both density grids
        pts = np.vstack([self.grids[n].grid_coords for n in range(len(
                                                                self.grids))])
        #and a concatenated list of the densities calculated for
        #both density grids
        vals = np.hstack([self.grids[n]._calc_density(x,
                                y).flatten() for n in range(len(self.grids))])

        #then interpolate from those points and values to the centerpoints
        #of all of the land centerpoints
        dens = interpolate.griddata(pts, vals, (self.land_gj, self.land_gi),
                                                            method = 'cubic')

        return dens


class _DirectionalitySurface:
    def __init__(self, dir_lyr, mixture, approximation_len = 7500,
                                                    vm_distr_kappa = 12):
        #dimensions
        self.dim = dir_lyr.dim
        #resolution (i.e. cell-size); defaults to 1
        self.res = dir_lyr.res

        #set the default values, in case they're accidentally fed in
        #as None values in the params
        if approximation_len is None:
            approximation_len = 7500
        if vm_distr_kappa is None:
            vm_distr_kappa = 12
            #NOTE: when I eventually have a mix = True/False argument
            #(to make shallow gradients work better),
            #should probably set this to a lower value if mix == False

        self.lyr_num = dir_lyr.idx
        self.approximation_len = approximation_len
        self.surf = _make_directionality_surface(dir_lyr.rast, 
            mixture = mixture, approximation_len = self.approximation_len,
            vm_distr_kappa = vm_distr_kappa)

        assert self.approximation_len == self.surf.shape[2], ("_"
            "DirectionalitySurface.approximation_len not equal to "
            "DirectionalitySurface.surf.shape[2]")

    def _draw_directions(self, x, y):
        choices = r.randint(low = 0, high = self.approximation_len,
                                                                size = len(x))
        return self.surf[y, x, choices]


class _KDTree:
    def __init__(self, coords, leafsize = 100):
        self.tree = cKDTree(data = coords, leafsize = leafsize)

    def _query(self, coords, dist, k=2):
        dists, pairs = self.tree.query(x = coords, k = k,
                                                    distance_upper_bound=dist)
        return(dists, pairs)


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#create strings from input cell coordinates
def _make_cell_strings(gi, gj, dim_om):
    #get strings for both i and j cooridnates, zfilling to the
    #correct order-of-magnitude (of the larger land dimensions)
    i_strs = [str(int(i)).zfill(dim_om) for i in gi.flatten()]
    j_strs = [str(int(j)).zfill(dim_om) for j in gj.flatten()]

    #join those strings to one single string, unique for each cell
    cells = [''.join(c) for c in list(zip(i_strs,j_strs))]

    return cells


#make a density grid, based on the Landscape object, the chosen window-width, 
#and the Boolean arguments dictating whether or not the grid's x- and y-
#dimension cells should be centered on the land edges (i.e. 0 and dim[_])
def _make_density_grid(land, ww, x_edge, y_edge):

    #half-window width
    hww = ww/2.

    #get land dimensions
    dim = land.dim
    dim_om = land._dim_om

    #create a dictionary of cell ranges, one for when cells center
    #on edge values (i.e. 0 and dim[n] for either dimension), 
    #the other for when they don't (i.e. run from hww to dim[n] - hww)
    edge_range_dict = {True:  np.arange(0, dim[0]+ww, ww),
                       False: np.arange(0+hww, dim[0]+hww, ww)}

    #create the meshgrid of the centerpoints of neighborhoods (or cells)
    #within which population will be counted
    #(x_edge and y_edge arguments determine whether this grid's 
    #x and y cells are centered on the landscape edges or not)
    #NOTE: these are expressed as points in continuous space from 0 to
    #each land dimension, NOT as cell
    #numbers (which will be calculated below)
    gj, gi = np.meshgrid(edge_range_dict[x_edge], edge_range_dict[y_edge])

    #and get flattened lists of the grid's i and j coordinates 
    j = gj.flatten()
    i = gi.flatten()


    #create a single, large Polygon object of the landscape quadrilateral
    land_poly_coords = ((0,0), (dim[0], 0), (dim[0], dim[1]), (0, dim[1]))
    land_poly = g.Polygon(land_poly_coords)


    #create a list of quadrilaterals centered on each of the points in the grid
    polys = [g.Polygon(((j[n]-hww, i[n]-hww), (j[n]-hww, i[n]+hww),
        (j[n]+hww, i[n]+hww), (j[n]+hww, i[n]-hww))) for n in range(len(j))]

    #use the Polygons' intersection method to calculate the total area
    #of the land that is covered by each grid cell (which will be 
    #used as the denominator for calculating neighborhood population densities
    #from neighborhood population counts)
    areas = np.reshape([p.intersection(land_poly).area for p in polys],
                                                                    gj.shape)

    #get lists of the integer-identifiers (i,j in the proper matrix sense)
    #of the cells in each of the grids (so that these can be matched up
    #with individuals' cell numbers when calculating population density while
    #the model is running)
    #(e.g. a cell centered at point y = 0, x = 5 would be cell i = 0,j = 6
    #if the window-width is 1, or cell i =0,j = 1 if the window width is 5)
    i_cells = (i - (hww * (y_edge))) // ww + (y_edge)
    j_cells = (j - (hww * (x_edge))) // ww + (x_edge)


    #turn the cell-number integers into cell strings
    #(NOTE: the previous algorithm used to calculate population density
    #was more or less the same, but found individuals' cells by
    #checking numerically whether they were within each window, but
    #was much slower and scaled poorly with decreasing window width,
    #increasing landscape size, and increasing population size.
    #This approach instead uses the floor-divide // on individuals'
    #x and y coordinates to generate cell-number integers for them,
    #converts those to strings, then counts the number of instances of each
    #string for the population and uses those values as the counts
    #of individuals within each density-grid cell,
    #obviating the need to loop across grid dimensions.
    #This performs better and scales much better.
    cells = _make_cell_strings(i_cells, j_cells, dim_om)

    #use the above-created data structures to create two DensityGrid objects
    #(which will inhere to the Landscape object as attributes)
    grid = _DensityGrid(dim, dim_om, ww, gi, gj, cells, areas,
                                            x_edge= x_edge, y_edge = y_edge)
    return grid


#create 4 density grids, one for each offset (i.e. each
#combination of offset by 0 and by 0.5*window_width)
def _make_density_grids(land, ww):
    #make a grid for each combo of Booleans for x_edge and y_edge
    g1 = _make_density_grid(land, ww, x_edge = True, y_edge = True)
    g2 = _make_density_grid(land, ww, x_edge = False, y_edge = False)
    g3 = _make_density_grid(land, ww, x_edge = True, y_edge = False)
    g4 = _make_density_grid(land, ww, x_edge = False, y_edge = True)
    return(g1, g2, g3, g4)


# Function to generate a simulative Von Mises mixture distribution
#sampler function
def _make_von_mises_mixture_sampler(neigh, dirs, vm_distr_kappa=12,
                                                    approximation_len = 5000):
    # Returns a lambda function that is a quick and reliable way to simulate
    #draws from a Von Mises mixture distribution:
    # 1.) Chooses a direction by neighborhood-weighted probability
    # 2.) Makes a random draw from a Von Mises dist centered on the direction,
    #     with a vm_distr_kappa value set such that the net effect, 
    #when doing this a ton of times for a given neighborhood and then
    #plotting the resulting histogram, gives the visually/heuristically
    #satisfying approximation of a Von Mises mixture distribution

    # NOTE: Just visually, heuristically, vm_distr_kappa = 10 seemed
    #like a perfect middle value (vm_distr_kappa ~3 gives too
    # wide of a Von Mises variance and just chooses values around
    #the entire circle regardless of neighborhood
    # probability, whereas vm_distr_kappa ~20 produces noticeable reductions
    #in probability of moving to directions between the 8 queen's
    #neighborhood directions (i.e. doesn't simulate the mixing well enough)
    #and would generate artefactual movement behavior); 12 also
    #seemed to really well in generating probability valleys
    # when tested on neighborhoods that should generate bimodal distributions
    d = list(dirs.ravel())
    n = list(neigh.ravel())
    del d[4]
    del n[4]
    sum_n = float(sum(n))
    if sum_n > 0:
        n_probs = [i / sum_n for i in n]
    else:
        n_probs = [.125]*8
    loc_choices = r.choice(d, approximation_len, replace = True, p = n_probs)
    loc_choices = list(C(loc_choices).items())
    approx = np.hstack([s_vonmises.rvs(vm_distr_kappa, loc=loc, scale=1,
                                size = size) for loc, size in loc_choices])
    return approx


# Runs the Von Mises mixture sampler function (_make_von_mises_mixture_sampler)
# or the Von Mises unimodal sampler function (make_von_mises_unimodal_sampler)
#across the entire landscape and returns an array-like (list of lists) of the 
#resulting lambda-function samplers
def _make_directionality_surface(rast, mixture=True, approximation_len=5000,
                                                        vm_distr_kappa=12):
    queen_dirs = np.array([[-3 * pi / 4, -pi / 2, -pi / 4], [pi, np.NaN, 0],
                                            [3 * pi / 4, pi / 2, pi / 4]])

    #copy the landscape raster
    rast = deepcopy(rast)

    # create embedded raster (so that the edge probabilities are
    #appropriately calculated)
    embedded_rast = np.zeros(shape=[n + 2 for n in rast.shape])
    embedded_rast[1:embedded_rast.shape[0
                                    ] - 1, 1:embedded_rast.shape[1] - 1] = rast

    #create a numpy array and store vectors approximating the functions!
    dir_surf = np.float16(np.zeros((rast.shape[0],
                                        rast.shape[1], approximation_len)))

    for i in range(rast.shape[0]):
        for j in range(rast.shape[1]):
            neigh = embedded_rast[i:i + 3, j:j + 3].copy()
            if mixture:
                dir_surf[i, j, :] = _make_von_mises_mixture_sampler(neigh,
                    queen_dirs, vm_distr_kappa = vm_distr_kappa,
                    approximation_len= approximation_len)
            else:
                dir_surf[i, j, :] = _make_von_mises_unimodal_sampler(neigh,
                    queen_dirs, vm_distr_kappa = vm_distr_kappa,
                    approximation_len= approximation_len)
    return dir_surf


#coarse wrapper around the nlmpy package
def _make_nlmpy_raster(nlmpy_params):
    if 'nlmpy' in globals():
        #pop out the name of the function to be called
        fn_name = nlmpy_params.pop('function')
        #try to create the nlmpy raster
        try:
            #get the function to be called
            fn = getattr(nlmpy, fn_name)
            nlm = fn(**nlmpy_params)
        except Exception as e:
            raise ValueError(("NLMpy could not generate the raster using the "
                "parameters provided. It threw the following error:\n\n\t"
                "%s\n\n.") % e)
        #if the nlm generated is not constrained between 0 and 1, rescale
        #to that range
        if nlm.min() < 0 or nlm.max() > 1:
            nlm, min_inval, max_inval = _scale_raster(nlm)
        return nlm
    else:
        raise ModuleNotFounderror(("It appears that 'nlmpy.nlmpy' was not "
            "imported successfully."))


#linearly scale a raster to 0 <= x <= 1, and return the min and max input vals
#as well (for possible reversion)
def _scale_raster(rast, min_inval=None, max_inval=None, min_outval=0,
                                                                max_outval=1):
    if min_inval is None:
        min_inval = rast.min()
    if max_inval is None:
        max_inval = rast.max()
    scale_rast = (rast - min_inval)/(max_inval - min_inval)
    return(scale_rast, min_inval, max_inval)

