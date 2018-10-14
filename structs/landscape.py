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
from ops import change

#other imports
import numpy as np
import numpy.random as r
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter as C
from copy import deepcopy
from operator import itemgetter as ig
from scipy import interpolate
from shapely import geometry as g


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

class Scape:

    #######################
    ### SPECIAL METHODS ###
    #######################

    def __init__(self, rast, scape_type, name, dim, res=(1,1), ulc = (0,0),
                 prj=None, scale_min=0, scale_max=1):
        self.idx = None
        self.type = scape_type
        self.name = str(name)
        self.dim = dim
        assert type(self.dim) in [tuple, list], "dim must be expressed on a tuple or a list"
        #set the resoultion (res; i.e. cell-size) and upper-left corner (ulc)
        #to defaults; will be reset if landscape read in from a GIS raster or
        #numpy txt array file
        self.res = res
        self.ulc = ulc
        self.prj = prj
        self.rast = rast
        assert type(self.rast) == np.ndarray, "rast should be a numpy.ndarray"
        self._scale_min = scale_min
        self._scale_max = scale_max


    #####################
    ### OTHER METHODS ###
    #####################

    #################
    #private methods#
    #################

    #method for writing the scape's raster to a geotiff raster file
    def _write_geotiff(self, filepath):
        io._write_geotiff(filepath, self)

    #method for writing the scape's raster to a numpy txt array
    def _write_txt_array(self, filepath):
        io._write_txt_array(filepath, self)


    ################
    #public methods#
    ################

    def plot(self, colorbar=True, im_interp_method='nearest', cmap = 'terrain', x=None, y=None, zoom_width=None, vmin=None, vmax=None):
        plt_lims = viz.get_plt_lims(self, x, y, zoom_width)
        viz.plot_rasters(self, colorbar = colorbar, im_interp_method = im_interp_method, cmap = cmap, plt_lims = plt_lims, vmin = vmin, vmax = vmax)

    #method for writing the scape's raster to a file of the specified format
    def write_raster(self, filepath, raster_format):
        assert raster_format in ['geotiff', 'txt'], ("The raster_format "
            "must be one of the following: 'geotiff', 'txt'.")
        if raster_format == 'geotiff':
            self._write_geotiff(filepath)
        elif raster_format == 'txt':
            self._write_txt_array(filepath)


class Land(dict):

    #######################
    ### SPECIAL METHODS ###
    #######################

    def __init__(self, scapes, res=(1,1), ulc=(0,0), prj=None):
        #check the scapes dict is correct, then update the Land with it
        assert False not in [scape.__class__.__name__ == 'Scape' for scape in scapes.values()], 'All landscapes supplied in scapes must be of type landscape.Scape.'
        self.update(scapes)
        #set the scapes attribute to the dict values (it still changes dynamically as things update, it runs
        #3x faster than calling values() during the model, and it makes code nicer looking
        #FIXME: I had been using the following lines, but deepcopy breaks because it cannot pickle dict.values() objects
        #self.scapes = self.values()
        #set the number of scapes in the stack
        self.n_scapes = len(self)
        #set the idx attributes on the scapes
        [setattr(v, 'idx', k) for k, v in self.items()]

        #check that scapes' dimensions are all the same
        assert len(set([land.dim for land in list(self.values())])) == 1, 'Dimensions of all landscapes must be equal.'
        #then set the dim to the 0th landscape's dim 
        self.dim = list(self.values())[0].dim
        #get the order of magnitude of the larger land dimension (used when zero-padding cell-coorindate strings)
        self._dim_om = max([len(str(d)) for d in self.dim])
        #set the resoultion (res; i.e. cell-size), upper-left corner (ulc), and
        #projection (prj) to the provided values
        self.res = res
        self.ulc = ulc
        self.prj = prj
        #TODO: DO I ACTUALLY HAVE TO RUN THIS LINE?
        #self.set_scapes_res_ulc_prj()
        #And check that the Scapes' res, ulc, and prj values are equal to 
        #the Land's values
        assert np.all([scape.res == self.res for scape in self.values()]), ("Not"
            " all Scapes have the same resolution value (attribute 'res') "
            "as the 'res' value being used to create the Land object that "
            "should contain them.")
        assert np.all([scape.ulc == self.ulc for scape in self.values()]), ("Not"
            " all Scapes have the same upper-left corner value (attribute "
            " 'ulc') as the 'ulc' value being used to create the Land object "
            "that should contain them.")
        assert np.all([scape.prj == self.prj for scape in self.values()]), ('Not'
            " all Scapes have the same projection (attribute 'prj') "
            "as the 'prj' value being used to create the Land object that "
            "should contain them.")

        #create a changer attribute (defaults to None, but will later be set
        #to an ops.change.LandChanger object if params call for it)
        self._changer = None

    #define the __str__ and __repr__ special methods
    #NOTE: this doesn't excellently fit the Python docs' specification for 
    #__repr__; I should massage this some more when done writing the codebase
    def __str__(self):
        #get the longest scape name, to be used to horizontally rectify 
        #all the lines for each of the scapes
        max_len_scape_name = max([len(scape.name) for scape in self.values()])
        #get a string representation of the class
        type_str = str(type(self))
        #get a string of all the scapes
        scapes_str = '%i Scape%s:\n' % (self.n_scapes, 's' * (len(self) > 1))
        scapes_str = scapes_str + '\n'.join(['\tscape %i: ' % k +
            "%s'%s':\t" % (' ' * (max_len_scape_name - len(v.name)), v.name) +
            str(v) for k,v in self.items()])
        #get a string representation of the first two and last two parameters
        params_str = "\nParameters:\n\t" + ',\n\t'.join(sorted([str(k) +
            ': ' + str(v) for k,v in vars(self).items()][:2])) + ','
        params_str = params_str + '\n\t...\n\t'
        params_str = params_str + ',\n\t'.join(sorted([str(k) +
            ': ' + str(v) for k,v in vars(self).items()][-2:]))

        return '\n'.join([type_str, scapes_str, params_str])

    def __repr__(self):
        repr_str = self.__str__()
        return repr_str


    #####################
    ### OTHER METHODS ###
    #####################

        #################
        #private methods#
        #################

    #method to set a raster
    def _set_raster(self, scape_num, rast):
        self[scape_num].rast = rast

    #method to res all scapes' res and ulc attributes 
    def _set_scapes_res_ulc_prj(self):
        [setattr(scape, 'res', self.res) for scape in self.values()]
        [setattr(scape, 'ulc', self.ulc) for scape in self.values()]
        [setattr(scape, 'prj', self.prj) for scape in self.values()]

    #method to make landscape changes
    def _make_change(self, t):
        self._changer._make_change(t)


        ################
        #public methods#
        ################

    #method to plot the landscape (or just a certain scape)
    def plot(self, scape_num=None, colorbar=True, cmap='terrain', im_interp_method='nearest', x=None, y=None, zoom_width=None, vmin=None, vmax=None):
        plt_lims = viz.get_plt_lims(self, x, y, zoom_width)
        viz.plot_rasters(self, scape_num = scape_num, colorbar = colorbar, im_interp_method = im_interp_method, cmap = cmap, plt_lims = plt_lims, vmin = vmin, vmax = vmax)

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

def _make_random_scape(dim, n_pts, interp_method="cubic", num_hab_types=2, dist='beta', alpha=0.05, beta=0.05):  # requires landscape to be square, such that dim = domain = range
    # NOTE: can use "nearest" interpolation to create random patches of habitat (by integer values); can change num_hab_types to > 2 to create a randomly multi-classed landscape
    # NOTE: I guess this could be used for rectangular landscapes, if the square raster is generated using the larger of the two dimensions, and then the resulting array is subsetted to the landscape's dimensions
    # NOTE: This seems to generate decent, believable random rasters!
    # n_pts/dim ratio:
    # ~0.01-0.05 --> broad, simple gradients
    # ~0.05-0.10 --> slightly more complex, occasionally landscape of low regions divided by a major high region (or vice versa)
    # ~0.10-0.50 --> landscape can be broken up into numerous fragmented patches (though sometimes still relatively homogeneous, with one or two small, extremely different patches
    # ~0.50-1.00 --> more highly fragmented version of above
    max_dim = max(dim)
    if interp_method == 'nearest':
        if dist == 'unif':
            vals = r.rand(n_pts) * (num_hab_types - 1)
        elif dist == 'beta':
            vals = r.beta(alpha, beta, n_pts) * (num_hab_types - 1)
    else:
        if dist == 'unif':
            vals = r.rand(n_pts)
        elif dist == 'beta':
            vals = r.beta(alpha, beta, n_pts)
    pts = r.normal(max_dim / 2, max_dim * 2, [n_pts,
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


def _make_defined_scape(dim, pts, vals, interp_method="cubic", num_hab_types=2):  # pts should be provided as n-by-2 Numpy array, vals as a 1-by-n Numpy array

    # NOTE: There seem to be some serious issues with all of this code, because the resulting landscapes are not quite symmetrical; and making small tweaks (e.g. adding 5 to all input points' coordinates) doesn't just tweak the output landscape but instead completely gets ride of it; I have an intuition that it comes from the code that coerces all raster values to 0 <= val <= 1, becuase that doesn't look it does quite what I originally intended for it to do, but I'm not really sure... anyhow, for now it works for my initial testing purposes

    # NOTE: like the _make_random_scape function, this also requires landscape to be square, but I guess this could be used for rectangular landscapes, if the square raster is generated using the larger of the two dimensions, and then the resulting array is subsetted to the landscape's dimensions

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


def _make_land(params, num_hab_types=2):
    # NOTE: If a multi-class (rather than binary) block-habitat raster would be of interest, would need to make
    # num_hab_types customizable)

    dim = params.land.main.dim
    res = params.land.main.res
    ulc = params.land.main.ulc
    prj = params.land.main.prj
    if res is None:
        res = (1,1)
    if ulc is None:
        ulc = (0,0)
    #leave default projection as None for now
    if prj is None:
        prj = None

    #create a dictionary to hold all the scapes to be created
    scapes = {}

    #create a list to hold the information for any rasters to be created from
    #file (because they will all be created after the loop over the scapes, so
    #that it's simpler to check agreement among raster resolutions and 
    #registrations
    file_scape_params = {'names': [], 'scape_nums':[], 'filepaths':[], 'scale_min_vals':[], 'scale_max_vals':[]}

    #then loop over the scapes in params.land.scapes and create each one
    for n, (scape_name, scape_params) in enumerate(params.land.scapes.items()):

        #get the init parameters
        init_params = deepcopy(scape_params.init)

        #determine which type of scape this is to be (valid: 'rand', 'defined', 'file')
        init_keys = [*init_params]
        if len(init_keys) > 1:
            raise ValueError(("The %ith scape (params['land']['scapes'] "
                "key '%s') appears to have parameters for more than one scape "
                "type.  Choose a single scape type (valid values: 'rand', "
                "'defined', 'file', 'nlmpy') and provide a sub-dictionary of parameters "
                "for only that type.") % (n, str(k)))
        scape_type = init_keys[0] 
        assert scape_type in ['rand', 'defined', 'file', 'nlmpy'], ("The "
            "parameters sub-dictionary for the %ith scape (params['land']"
            "['scapes'] key '%s') has an invalid key value. Valid keys are: "
            "'rand', 'defined', 'file'.") % (n, str(k))

        #create a random scape, if called for
        if scape_type == 'rand':
            #create the random scape and add it to the scapes dict
            scape_rast = _make_random_scape(dim, **init_params[scape_type], 
                                           num_hab_types=num_hab_types)
            scapes[n] = Scape(scape_rast, scape_type = scape_type, 
                name = scape_name, dim = dim, res = res, ulc = ulc)

        #or else create a defined scape 
        elif scape_type == 'defined':
            #create the defined raster
            scape_rast = _make_defined_scape(dim, **init_params[scape_type],
                                            num_hab_types=num_hab_types)
            scapes[n] = Scape(scape_rast, scape_type = scape_type, 
                name = scape_name, dim = dim, res = res, ulc = ulc)

        #or else create an nlmpy scape
        elif scape_type == 'nlmpy':
            #get the params
            nlmpy_params = init_params[scape_type]
            #make the nlm
            scape_rast = spt._make_nlmpy_raster(nlmpy_params)
            #check that its dimensions match those of the Land
            assert scape_rast.shape == dim, ("The dimensions of the NLM created"
                " appear to differ from the Land dims: Land has dims %s, NLM "
                "has dims %s.\n\n") % (str(dim), str(scape_rast.shape))
            #set the nlm as the nth Scape in the Land
            scapes[n] = Scape(scape_rast, scape_type = scape_type,
                name = scape_name, dim = dim, res = res, ulc = ulc)

        #or else create a scape from file
        elif scape_type == 'file':
            #set this scape to None, temporarily
            scapes[n] = None
            #then store the relevant info for this scape in the file_scape_params dict; this scape's raster
            #will be replaced with a GIS or numpy raster after the Land is #created (NOTE: easier to do this all at 
            #once afterward, to check that the resolution and registration of the rasters all agree)
            file_scape_params['scape_nums'].append(n)
            file_scape_params['names'].append(scape_name)
            [file_scape_params[k+'s'].append(v) for k,v in init_params[scape_type].items()]

    #now set the necessary layers to their file rasters, if applicable
    if True in [len(v) > 0 for v in file_scape_params.values()]:
        file_scapes, res, ulc, prj = get_file_rasters(land_dim = dim, **file_scape_params)
        for n in file_scape_params['scape_nums']:
            scapes[n] = file_scapes[n]

    #create the land object
    land = Land(scapes, res=res, ulc=ulc, prj=prj)

    #grab the change parameters into a dictionary of scape_num:events:events_params hierarchical
    change_params = {k:v.change for k,v in params.land.scapes.items() if 'change' in v.keys()} 
    #and then replace the land.changer attribute with an ops.change.LandChanger object, if params call for it
    if len(change_params) > 0:
        land._changer = change._LandChanger(land, change_params)

    return land


def _get_file_rasters(land_dim, names, scape_nums, filepaths, scale_min_vals, scale_max_vals):
    assert len(scape_nums)==len(filepaths), 'Parameters provide a different number of GIS raster files to read in than of scape numbers to set them to .'
    res = []
    ulc = []
    prj = []
    rasters = []
    scapes = []
    for n,filepath in enumerate(filepaths):
        #get the array, dim, res, upper-left corner, and projection from io.read_raster
        rast_array, rast_dim, rast_res, rast_ulc, rast_prj = io._read_raster(filepath, dim)
        #check that the dimensions are right
        assert rast_dim == land_dim, 'Variable land_dim and the dimensions of the input raster %s appear to differ. Please clip %s to the correct dimensions and try again. %s has dimensions (%i, %i), but land has dimensions (%i,%i)' % (filepath, filepath, filepath, rast_dim[0], rast_dim[1], land_dim[0], land_dim[1])
        #scale the raster to 0<=x<=1
        rast_array, scale_min, scale_max = spt._scale_raster(rast_array, min_inval = scale_min_vals[n], max_inval = scale_max_vals[n])
        #add the rast_array to the rast_arrays list
        rasters.append(rast_array)
        #and add the raster's ulc and res attributes to the appropriate lists
        res.append(tuple(rast_res))
        ulc.append(tuple(rast_ulc))
        prj.append(rast_prj)
        #and update scale_min_vals and scale_max_vals if the values used in spt.scale_raster() don't equal the
        #values from params (should only occur if the params values provided were None)
        if scale_min != scale_min_vals[n]:
            scale_min_vals[n] = scale_min
        if scale_max != scale_max_vals[n]:
            scale_max_vals[n] = scale_max
    #check that there is only one unique tuple each in res and ulc lists,
    #and only one unique WKT string in the prj list
    #(otherwise rasters are not in the same projection and/or registration)
    res_C = C(res)
    ulc_C = C(ulc)
    prj_C = C(prj)
    assert len(res_C) == 1 and list(res_C.values())[0] == len(res), ('Some of '
            'the rasters read in from file appear to differ in resolution.')
    assert len(ulc_C) == 1 and list(ulc_C.values())[0] == len(ulc), ('Some of '
            'the rasters read in from file appear to differ in registration (i.e. '
            'top-left corner).')
    assert len(prj_C) == 1 and list(prj_C.values())[0] == len(prj), ('Some of '
            'the rasters read in from file appear to have different projections.')
    #then get the res and ulc values
    res = res[0]
    ulc = ulc[0]
    prj = prj[0]
    #create scapes from the rasters
    scapes = [Scape(rast, scape_type = 'file', name = name, dim = land_dim,
            res = res, ulc = ulc, prj = prj, scale_min=scale_min,
            scale_max=scale_max) for scape_num, name, rast, scale_min,
            scale_max in zip(scape_nums, names, rasters, scale_min_vals, scale_max_vals)]
    return(scapes, res, ulc, prj)


# function for reading in a pickled landscape stack
def read_pickled_land(filename):
    import cPickle
    with open(filename, 'rb') as f:
        land = cPickle.load(f)
    return land

