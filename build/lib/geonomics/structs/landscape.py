#!/usr/bin/python
# landscape.py

'''
Defines the Layer and Landscape classes, with their associated methods and
supporting functions
'''
#geonomics imports
from geonomics.utils.viz import _get_plt_lims, _plot_rasters, _check_display
from geonomics.utils.io import _write_geotiff, _write_txt_array, _read_raster
from geonomics.utils.spatial import _scale_raster, _make_nlmpy_raster
from geonomics.ops.change import _LandscapeChanger

#other imports
import numpy as np
import numpy.random as r
import random
import matplotlib as mpl
_check_display()
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

class Layer:
    """
    Representation of a single environmental layer (i.e. variable).

    Multiple Layers are collected as serial integer-keyed values within
    a Landscape dict.

    Attributes
    ----------

        NOTE: For more detail, see the documentation for the parameters
              that correspond to many of the following attributes.

        coord_prec:
            The precision (i.e. number of significant digits) to which
            coordinate values should be round when the Layer is plotted.
       
        dim:
            The x,y (i.e. lon,lat; or i,j in array terms) dimensions of the
            Layer. (Must be the same as the Landscape to which the Layer
            belongs.)

        idx:
            Index number of the Layer (i.e. its key within the Landscape dict)

        name:
            The string name of the Layer

        prj:
            The projection of the Layer (formatted as a PROJ4 string).
            (Must be the same as the Landscape to which the Layer belongs.)

        rast:
            The 2d numpy array, of shape `Landscape.dim`, containing the
            environmental values for this Layer.

        res:
            The x,y (i.e. lon,lat; or i,j in array terms) spatial resolution
            (i.e. cell sizes) of the layer.  (Must be the same as the
            Landscape to which the Layer belongs.)


        type:
            A string indicating the type of the Layer ('random', 'defined',
            'file', or 'nlmpy')

        ulc: 
            The x,y (i.e. lon,lat; or i,j in array terms) coordinates of the
            upper left corner of the layer. (Must be the same as the Landscape
            to which the Layer belongs.)

        units:
            A string representation of the units of the Layer's variable
            (to be used for plotting).
    
    """

    #######################
    ### SPECIAL METHODS ###
    #######################

    def __init__(self, rast, lyr_type, name, dim, res=(1,1), ulc=(0,0),
                 prj=None, coord_prec=0, units='', scale_min=0, scale_max=1):
        self.idx = None
        self.type = lyr_type
        self.name = str(name)
        #the x,y (i.e. j,i) dimensions
        self.dim = dim
        assert type(self.dim) in [tuple, list], ("dim must be expressed "
                                                 "as a tuple or a list")
        #set the resoultion (res; i.e. x,y cell sizes) and upper-left corner
        #(ulc) to defaults; will be reset if layer read in from a GIS raster or
        #numpy txt array file
        self.res = res
        self.ulc = ulc
        self.prj = prj
        self.coord_prec = coord_prec
        self.units = units
        # attributes for arrays of cell-bounds coordinates, to be set when a
        # Landscape is instantiated
        self._x_cell_bds = None
        self._y_cell_bds = None
        self.rast = rast
        assert type(self.rast) == np.ndarray, "rast should be a numpy.ndarray"
        self._scale_min = scale_min
        self._scale_max = scale_max
        #attribute that will track the Species for which this Layer
        #is the K_layer (to be used to update
        #Species.K if this Layer undergoes any landscape changes)
        self._is_K = []


    #####################
    ### OTHER METHODS ###
    #####################

    #################
    #private methods#
    #################

    #method for writing the layer's raster to a geotiff raster file
    def _write_geotiff(self, filepath):
        _write_geotiff(filepath, self)

    #method for writing the layer's raster to a numpy txt array
    def _write_txt_array(self, filepath):
        _write_txt_array(filepath, self)

    # method for the getting the vmin and vmax values for plotting
    def _get_plot_vmin_vmax(self):
        return (self._scale_min, self._scale_max)

    #method for plotting the layer
    def _plot(self, cbar=True, cmap=None, x=None, y=None, zoom_width=None,
              vmin=None, vmax=None, ticks=None, mask_rast=None):
        plt_lims = _get_plt_lims(self, x, y, zoom_width)
        # get the vmin and vmax values, if not provided
        if vmin is None and vmax is None:
            vmin, vmax = self._get_plot_vmin_vmax()
        _plot_rasters(self, cbar=cbar, cmap=cmap,
                      plt_lims=plt_lims, vmin=vmin, vmax=vmax, ticks=ticks,
                      mask_rast=mask_rast)

    # method to get the tickmarks to use to plot the layer
    def _get_coord_ticks(self):
        x_ticks = np.int64(np.round(np.linspace(0, self.dim[0]-1, 8)))
        y_ticks = np.int64(np.round(np.linspace(0, self.dim[1]-1, 8)))
        x_tick_labs = [round(self.ulc[0] + (self.res[0] * (loc + 0.5)),
                             min(3, self.coord_prec)) for loc in x_ticks]
        y_tick_labs = [round(self.ulc[1] + (self.res[1] * (loc + 0.5)),
                             min(3, self.coord_prec)) for loc in y_ticks]
        return x_ticks, x_tick_labs, y_ticks, y_tick_labs

    # method for recovering a file-based raster's values expressed in its
    # native units (i.e. undoing the 0-1 scaling)
    def _get_rast_in_native_units(self):
        native = (self.rast * (self._scale_max
                               - self._scale_min)) + self._scale_min
        return native

    # method for getting colorbar ticks
    def _get_cbar_ticks_and_minmax_scaled_vals(self):
        #if self.type == 'file':
            #rast = self._get_rast_in_native_units()
            #min_val = rast.min()
            #max_val = rast.max()
        #else:
        #    min_val = self.rast.min()
        #    max_val = self.rast.max()
        min_val = self._scale_min
        max_val = self._scale_max
        ticks = np.round(np.linspace(min_val, max_val, 5), 3)
        return ticks, min_val, max_val


    #method for writing the lyr's raster to a file of the specified format
    def _write_raster(self, filepath, raster_format):
        assert raster_format in ['geotiff', 'txt'], ("The raster_format "
            "must be one of the following: 'geotiff', 'txt'.")
        if raster_format == 'geotiff':
            self._write_geotiff(filepath)
        elif raster_format == 'txt':
            self._write_txt_array(filepath)


class Landscape(dict):
    """
    Representation of a multi-layer (i.e. multivariate) landscape.

    Organized as a dict of multiple, serial integer-keyed Layer objects,
    
    Because the Landscape class inherits from `dict`, Layers can
    be indexed out using their index-number keys (e.g. `mod.land[<idx>]`).

    The Landscape is stored as the 'mod.land' attribute of its corresponding
    Model object.

    Attributes
    ----------

        NOTE: For more detail, see the documentation for the parameters
              that correspond to many of the following attributes.

        dim:
            The x,y (i.e. lon,lat; or i,j in array terms) dimensions of the
            Layer. (Must be the same as the Landscape to which the Layer
            belongs.)

        n_lyrs:
            The number of Layers in the Landscape

        prj:
            The projection of the Layer (formatted as a PROJ4 string).
            (Must be the same as the Landscape to which the Layer belongs.)

        res:
            The x,y (i.e. lon,lat; or i,j in array terms) spatial resolution
            (i.e. cell sizes) of the layer.  (Must be the same as the
            Landscape to which the Layer belongs.)

        ulc: 
            The x,y (i.e. lon,lat; or i,j in array terms) coordinates of the
            upper left corner of the layer. (Must be the same as the Landscape
            to which the Layer belongs.)
            
    """

    #######################
    ### SPECIAL METHODS ###
    #######################

    def __init__(self, lyrs, res=(1,1), ulc=(0,0), prj=None, mod=None):
        #check the lyrs dict is correct, then update the Landscape with it
        assert False not in [
            lyr.__class__.__name__ == 'Layer' for lyr in lyrs.values(
            )], ('All layers supplied in lyrs must be of type '
            'landscape.Layer.')
        self.update(lyrs)
        #set the lyrs attribute to the dict values (it still changes
        #dynamically as things update, it runs
        #3x faster than calling values() during the model, and it makes
        #code nicer looking
        #FIXME: I had been using the following lines, but deepcopy breaks
        #because it cannot pickle dict.values() objects
        #self.lyrs = self.values()
        #set the number of lyrs in the stack
        self.n_lyrs = len(self)
        #set the idx attributes on the lyrs
        [setattr(v, 'idx', k) for k, v in self.items()]

        #check that lyrs' dimensions are all the same
        assert len(set([land.dim for land in list(self.values(
            ))])) == 1, 'Dimensions of all layers must be equal.'
        #then set the dim to the 0th layer's dim 
        self.dim = list(self.values())[0].dim
        #get the order of magnitude of the larger land dimension
        #(used when zero-padding cell-coorindate strings)
        self._dim_om = max([len(str(d)) for d in self.dim])
        #set the resoultion (res; i.e. x,y cell sizes), upper-left corner
        #(ulc), and projection (prj) to the provided values
        self.res = res
        # get the resolution ratio (for calculating individuals' movement
        # distances on non-square-resolution rasters)
        self._res_ratio = tuple(np.abs([val/max(self.res) for val in
                                        self.res]))
        self.ulc = ulc
        self.prj = prj
        # get arrays of grid-cell x and y bounds , for plotting
        # with plt.pcolormesh
        x_cell_bds, y_cell_bds = [np.linspace(self.ulc[i],
                                              self.ulc[i] + (self.res[i] * (
                                                             self.dim[i])),
                                              self.dim[i] + 1) for i in range(
                                                                            2)]
        self._x_cell_bds = x_cell_bds
        self._y_cell_bds = y_cell_bds
        #And check that the Layers' res, ulc, and prj values are equal to 
        #the Landscape's values
        assert np.all([lyr.res == self.res for lyr in self.values(
            )]), ("Not all Layers have the same resolution value "
            "(attribute 'res') as the 'res' value being used to create "
            "the Landscape object that should contain them.")
        assert np.all([lyr.ulc == self.ulc for lyr in self.values(
            )]), ("Not all Layers have the same upper-left corner value "
            "(attribute 'ulc') as the 'ulc' value being used to create "
            "the Landscape object that should contain them.")
        assert np.all([lyr.prj == self.prj for lyr in self.values(
            )]), ("Not all Layers have the same projection (attribute 'prj') "
            "as the 'prj' value being used to create the Landscape "
            "object that should contain them.")
        #check that all layers' values are constrained to 0 <= x <= 1
        for lyr in lyrs.values():
            assert np.all(0 <= lyr.rast), ("Layer '%s' contains values "
                                           "less than 0.") %(lyr.name)
            assert np.all(lyr.rast <= 1), ("Layer '%s' contains values "
                                           "greater than 1.") %(lyr.name)


        #create a changer attribute (defaults to None, but will later be set
        #to an ops.change.LandscapeChanger object if params call for it)
        self._changer = None

    #define the __str__ and __repr__ special methods
    #NOTE: this doesn't excellently fit the Python docs' specification for 
    #__repr__; I should massage this some more when done writing the codebase
    def __str__(self):
        #get the longest lyr name, to be used to horizontally rectify 
        #all the lines for each of the lyrs
        max_len_lyr_name = max([len(lyr.name) for lyr in self.values()])
        #get a string representation of the class
        type_str = str(type(self))
        #get a string of all the lyrs
        lyrs_str = '%i Layer%s:\n' % (self.n_lyrs, 's' * (len(self) > 1))
        lyrs_str = lyrs_str + '\n'.join(['\tlyr %i: ' % k +
            "%s'%s':\t" % (' ' * (max_len_lyr_name - len(v.name)), v.name) +
            str(v) for k,v in self.items()])
        #get a string representation of the first two and last two parameters
        #params_str = "\nParameters:\n\t" + ',\n\t'.join(sorted([str(k) +
        #    ': ' + str(v) for k,v in vars(self).items()][:2])) + ','
        #params_str = params_str + '\n\t...\n\t'
        #params_str = params_str + ',\n\t'.join(sorted([str(k) +
        #    ': ' + str(v) for k,v in vars(self).items()][-2:]))
        #return '\n'.join([type_str, lyrs_str, params_str])
        return '\n'.join([type_str, lyrs_str])

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
    def _set_raster(self, lyr_num, rast):
        self[lyr_num].rast = rast

    #method to set all lyrs' res and ulc attributes
    def _set_lyrs_res_ulc_prj(self):
        [setattr(lyr, 'res', self.res) for lyr in self.values()]
        [setattr(lyr, 'ulc', self.ulc) for lyr in self.values()]
        [setattr(lyr, 'prj', self.prj) for lyr in self.values()]

    #method to make landscape changes
    def _make_change(self, t, verbose=False):
        self._changer._make_change(t=t, additional_args={'land': self},
                                   verbose=verbose)

    #method to plot the landscape (or just a certain lyr)
    def _plot(self, lyr_num=None, cbar=True, cmap=None, x=None, y=None,
              zoom_width=None, vmin=None, vmax=None, ticks=None,
              mask_rast=None):
        plt_lims = _get_plt_lims(self, x, y, zoom_width)
        # get the vmin and vmax values, if not provided
        if vmin is None and vmax is None:
            if lyr_num is None:
                vmins_vmaxs = [lyr._get_plot_vmin_vmax(
                                                    ) for lyr in self.values()]
                vmin = [v[0] for v in vmins_vmaxs]
                vmax = [v[1] for v in vmins_vmaxs]
            else:
                vmin, vmax = self[lyr_num]._get_plot_vmin_vmax()
        _plot_rasters(self, lyr_num=lyr_num, cbar=cbar, cmap=cmap,
                      plt_lims=plt_lims, vmin=vmin, vmax=vmax, ticks=ticks,
                      mask_rast=mask_rast)

    # method to get the tickmarks to use to plot the layer
    def _get_coord_ticks(self):
        try:
            coord_prec = min([lyr.coord_prec for lyr in self.values(
                            ) if lyr.coord_prec != 0])
        except ValueError:
            coord_prec = 0
        x_ticks = np.int64(np.round(np.linspace(0, self.dim[0]-1, 8)))
        y_ticks = np.int64(np.round(np.linspace(0, self.dim[1]-1, 8)))
        x_tick_labs = [round(self.ulc[0] + (self.res[0] * (loc + 0.5)),
                             min(3, coord_prec)) for loc in x_ticks]
        y_tick_labs = [round(self.ulc[1] + (self.res[1] * (loc + 0.5)),
                             min(3, coord_prec)) for loc in y_ticks]
        return x_ticks, x_tick_labs, y_ticks, y_tick_labs

        ################
        #public methods#
        ################

    # method for pickling a landscape
    def _write_pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

def _make_random_lyr(dim, n_pts, interp_method="cubic", num_hab_types=2,
                       dist='beta', alpha=0.05, beta=0.05):
    # requires layer to be square, such that dim = domain = range

    # NOTE: can use "nearest" interpolation to create random patches of
    #habitat (by integer values); can change num_hab_types to > 2
    #to create a randomly multi-classed layer
    # NOTE: I guess this could be used for rectangular layers, if the
    #square raster is generated using the larger of the two dimensions,
    #and then the resulting array is subsetted to the layer's dimensions
    # n_pts/dim ratio:
    # ~0.01-0.05 --> broad, simple gradients
    # ~0.05-0.10 --> slightly more complex, occasionally layer of
    #low regions divided by a major high region (or vice versa)
    # ~0.10-0.50 --> layer can be broken up into numerous fragmented
    #patches (though sometimes still relatively homogeneous, with
    #one or two small, extremely different patches
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
    # select seed points from well outside the eventaul layer, to ensure
    #interpolation across area of interest
    pts = r.normal(max_dim / 2, max_dim * 2, [n_pts,
                                              2])
    grid_x, grid_y = np.mgrid[1:max_dim:complex("%ij" % max_dim),
                            1:max_dim:complex( "%ij" % max_dim)]
        #NOTE: by this function's definition, the use of complex numbers
        #in here specifies the number of steps desired
    I = interpolate.griddata(pts, vals, (grid_x, grid_y), method=interp_method)
    if interp_method == 'nearest':
    # i.e., if being used to generate random habitat patches...
        I = I.round().astype(float)
    if interp_method == 'cubic':
        # transform to constrain all values to 0 <= val <= 1
        I = I + abs(I.min()) + (
        # NOTE: adding a small jitter to keep values from reaching == 0
        #or == 1, as would likely be the case with linear interpolation
                0.01 * r.rand())
        I = I / (I.max() + (0.01 * r.rand()))
    #use the dim tuple to subset an approriate size if dim not equal
    if dim[0] != dim[1]:
        I = I[:dim[0], :dim[1]]
    return I


def _make_defined_lyr(dim, rast, pts, vals, interp_method="cubic",
                                                    num_hab_types=2):
    #pts should be provided as n-by-2 Numpy array, vals as a 1-by-n Numpy array

    # NOTE: There seem to be some serious issues with all of this code,
    #because the resulting layers are not quite symmetrical;
    #and making small tweaks (e.g. adding 5 to all input points'
    #coordinates) doesn't just tweak the output layer but instead
    #completely gets ride of it; I have an intuition that it comes
    #from the code that coerces all raster values to 0 <= val <= 1,
    #becuase that doesn't look it does quite what I originally intended for
    #it to do, but I'm not really sure... anyhow, for now it works for
    #my initial testing purposes

    # NOTE: like the _make_random_lyr function, this also 
    #requires layer to be square, but I guess this could be
    #used for rectangular layers, if the square raster is
    #generated using the larger of the two dimensions, and
    #then the resulting array is subsetted to the layer's dimensions

    # NOTE: if discrete habitat patches desired, values should still
    #be fed in as proportions, and will then be multipled by
    #num_hab_types to develop habitat class values
    if rast is not None:
        I = rast
    else:
        if interp_method == 'nearest':
            vals = vals * (num_hab_types - 1)

        max_dim = max(dim)
        grid_x, grid_y = np.mgrid[1:max_dim:complex("%ij" % max_dim),
                                  1:max_dim:complex("%ij" % max_dim)]
        I = interpolate.griddata(pts, vals, (grid_x, grid_y), method=interp_method)
        if interp_method == 'nearest':
        # i.e., if being used to generate random habitat patches...
            I = I.round().astype(float)
        if interp_method == 'cubic':
            # transform to constrain all values to 0 <= val <= 1
            I = I + abs(I.min()) + (
                    0.01 * r.rand())
                # NOTE: adding a small jitter to keep values from
                #reaching == 0 or == 1,
            # as would likely be the case with linear interpolation
            I = I / (I.max() + (0.01 * r.rand()))
        #use the dim tuple to subset an approriate size if dim not equal
        if dim[0] != dim[1]:
            I = I[:dim[0], :dim[1]]
    return I


def _make_landscape(mod, params, num_hab_types=2, verbose=False):
    #print verbose output
    if verbose:
        print('\tMAKING LANDSCAPE...\n')
    # NOTE: If a multi-class (rather than binary) block-habitat raster
    #would be of interest, would need to make
    # num_hab_types customizable)
    dim = params.landscape.main.dim
    res = params.landscape.main.res
    ulc = params.landscape.main.ulc
    prj = params.landscape.main.prj
    if res is None:
        res = (1,1)
    if ulc is None:
        ulc = (0,0)
    #leave default projection as None for now
    if prj is None:
        prj = None

    #create a dictionary to hold all the lyrs to be created
    lyrs = {}

    #create a list to hold the information for any rasters to be created from
    #file (because they will all be created after the loop over the lyrs, so
    #that it's simpler to check agreement among raster resolutions and 
    #registrations
    file_lyr_params = {'names': [],
                       'lyr_nums': [],
                       'filepaths': [],
                       'scale_min_vals': [],
                       'scale_max_vals': [],
                       'coord_precs': [],
                       'unitss': []}

    #then loop over the lyrs in params.landscape.lyrs and create each one
    for n, (lyr_name, lyr_params) in enumerate(
                                        params.landscape.layers.items()):

        #get the init parameters
        init_params = deepcopy(lyr_params.init)

        #determine which type of lyr this is to be
        #(valid: 'random', 'defined', 'file')
        init_keys = [*init_params]
        if len(init_keys) > 1:
            raise ValueError(("The %ith layer (params['land']['layers'] "
                "key '%s') appears to have parameters for more than one layer "
                "type.  Choose a single layer type (valid values: 'random', "
                "'defined', 'file', 'nlmpy') and provide a sub-dictionary "
                "of parameters for only that type.") % (n, str(k)))
        lyr_type = init_keys[0]
        assert lyr_type in ['random', 'defined', 'file', 'nlmpy'], ("The "
            "parameters sub-dictionary for the %ith layer (params['land']"
            "['layers'] key '%s') has an invalid key value. Valid keys are: "
            "'random', 'defined', 'file'.") % (n, str(k))

        #create a random lyr, if called for
        if lyr_type == 'random':
            #create the random lyr and add it to the lyrs dict
            lyr_rast = _make_random_lyr(dim, **init_params[lyr_type],
                                           num_hab_types=num_hab_types)
            lyrs[n] = Layer(lyr_rast, lyr_type = lyr_type,
                name = lyr_name, dim = dim, res = res, ulc = ulc)

        #or else create a defined lyr 
        elif lyr_type == 'defined':
            #create the defined raster
            lyr_rast = _make_defined_lyr(dim, **init_params[lyr_type],
                                            num_hab_types=num_hab_types)
            lyrs[n] = Layer(lyr_rast, lyr_type = lyr_type,
                name = lyr_name, dim = dim, res = res, ulc = ulc)

        #or else create an nlmpy lyr
        elif lyr_type == 'nlmpy':
            #get the params
            nlmpy_params = init_params[lyr_type]
            #make the nlm
            lyr_rast = _make_nlmpy_raster(nlmpy_params)
            #check that its dimensions match those of the Landscape
            assert lyr_rast.shape == dim, ("The dimensions of the NLM "
                "created appear to differ from the Landscape dims: "
                "Landscape has x,y (i.e. j,i) dims %s, "
                "NLM has dims %s.\n\n") % (str(dim), str(lyr_rast.shape))
            #set the nlm as the nth Layer in the Landscape
            lyrs[n] = Layer(lyr_rast, lyr_type = lyr_type,
                name = lyr_name, dim = dim, res = res, ulc = ulc)

        #or else create a lyr from file
        elif lyr_type == 'file':
            #set this lyr to None, temporarily
            lyrs[n] = None
            #then store the relevant info for this lyr in the 
            #file_lyr_params dict; this lyr's raster
            #will be replaced with a GIS or numpy raster after the Landscape is 
            #created (NOTE: easier to do this all at 
            #once afterward, to check that the resolution and registration
            #of the rasters all agree)
            file_lyr_params['lyr_nums'].append(n)
            file_lyr_params['names'].append(lyr_name)
            [file_lyr_params[k+'s'].append(v) for k,v in init_params[
                                                        lyr_type].items()]

    #now set the necessary layers to their file rasters, if applicable
    if True in [len(v) > 0 for v in file_lyr_params.values()]:
        file_lyrs, res, ulc, prj, = _get_file_rasters(land_dim=dim,
                                                      **file_lyr_params)
        for n in file_lyr_params['lyr_nums']:
            lyrs[n] = file_lyrs[n]

        #set all other layers in the lyrs dict to have the same res, ulc,
        #and prj (because they all were not read in from file, whereas all
        #file-based layers were, and their res, ulc, and prj values were
        #already checked for equality, so we make the strong assumption
        #that the user wants the remaining non-file layers to be represented
        #by the same res, ulc, and prj values)
        for lyr_name in [*lyrs]:
            lyrs[lyr_name].res = res
            lyrs[lyr_name].ulc = ulc
            lyrs[lyr_name].prj = prj
        [setattr(lyrs[n], 'coord_prec',
                 v) for n, v in zip(file_lyr_params['lyr_nums'],
                                    file_lyr_params['coord_precs'])]

    #create the land object
    land = Landscape(lyrs, res=res, ulc=ulc, prj=prj, mod=mod)

    # set all Layers' cell-bounds attributes
    [setattr(lyr, '_x_cell_bds', land._x_cell_bds) for lyr in land.values()]
    [setattr(lyr, '_y_cell_bds', land._y_cell_bds) for lyr in land.values()]

    #grab the change parameters into a dictionary of
    #lyr_num:events:events_params hierarchical
    change_params = {k:v.change for k,v in params.landscape.layers.items(
                                                ) if 'change' in v.keys()}
    #and then replace the land._changer attribute with
    #an ops.change.LandscapeChanger object, if params call for it
    if len(change_params) > 0:
        #replace the keys, which are layer names, with layer numbers
        lyr_num_change_params = {}
        for k,v in change_params.items():
            lyr_num = [num for num, lyr in land.items() if lyr.name == k]
            assert len(lyr_num) == 1, ("Expected to find a single Layer "
                "matching the layer names for each Landscape-change event, "
                "but found %i Layers matching Layer name "
                "'%s'.") % (len(lyr_num), k)
            lyr_num_change_params[lyr_num[0]] = v
        land._changer = _LandscapeChanger(land, lyr_num_change_params,
            mod = mod)

    return land


def _get_file_rasters(land_dim, names, lyr_nums, filepaths, coord_precs,
                      scale_min_vals, scale_max_vals, unitss):
    assert len(lyr_nums)==len(filepaths), ('Parameters provide a '
                                           'different number of GIS raster '
                                           'files to read in than of layer '
                                           'numbers to set them to.')
    res = []
    ulc = []
    prj = []
    rasters = []
    lyrs = []
    for n,filepath in enumerate(filepaths):
        #get array, dim, res, ulc, and prj from io.read_raster
        rast_array, rast_dim, rast_res, rast_ulc, rast_prj = _read_raster(
            filepath, coord_precs[n], land_dim)
        #check that the dimensions are right
        assert rast_dim == land_dim, ('Variable land_dim and the dimensions '
            'of the input raster %s appear to differ. Please clip %s to '
            'the correct dimensions and try again, because it has x,y '
            'dimensions (%i, %i), but land has x,y '
            'dimensions (%i,%i)') % (filepath,
            filepath, rast_dim[0], rast_dim[1], land_dim[0], land_dim[1])
        #if either the scale_min_val or scale_max_val for this Layer is None,
        #set it to the min or max value of this Layer's array
        if scale_min_vals[n] is None:
            scale_min_vals[n] == rast_array.min()
        if scale_max_vals[n] is None:
            scale_max_vals[n] == rast_array.max()
        #scale the raster to 0<=x<=1
        rast_array, scale_min, scale_max = _scale_raster(rast_array,
                min_inval = scale_min_vals[n], max_inval = scale_max_vals[n])
        #add the rast_array to the rast_arrays list
        rasters.append(rast_array)
        #and add the raster's ulc and res attributes to the appropriate lists
        res.append(tuple(rast_res))
        ulc.append(tuple(rast_ulc))
        prj.append(rast_prj)
        #and update scale_min_vals and scale_max_vals if the values used
        #in spt.scale_raster() don't equal the values from params
        #(should only occur if the params values provided were None)
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
    #create lyrs from the rasters
    lyrs = [Layer(rast, lyr_type='file', name=name, dim=land_dim, res=res,
                  ulc=ulc, prj=prj, scale_min=scale_min, scale_max=scale_max,
                 units=units) for lyr_num, name,
            rast, scale_min, scale_max, units in zip(lyr_nums, names, rasters,
            scale_min_vals, scale_max_vals, unitss)]
    return(lyrs, res, ulc, prj)


# function for reading in a pickled landscape stack
def read_pickled_land(filename):
    import cPickle
    with open(filename, 'rb') as f:
        land = cPickle.load(f)
    return land

