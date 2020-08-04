#!/usr/bin/python
# change.py

'''
Classes and functions to implement landscape, demographic, and life-history
parameter change operations.
'''
#geonomics imports
from geonomics.utils import spatial as spt
from geonomics.utils import io
from geonomics.utils.viz import _check_display

#other imports
import os
import numpy as np
import numpy.random as r
import matplotlib as mpl
_check_display()
import matplotlib.pyplot as plt
from operator import attrgetter as ag
from collections import OrderedDict as OD
from collections import Counter as C
from copy import deepcopy


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

#base class for _LandscapeChanger and _SpeciesChanger
class _Changer:
    def __init__(self, params):
        #set the type-label of the changer
        self.type = None

        #grab a deep copy of the change params
        self.change_params = deepcopy(params)

        #create self.changes, which defaults to None but will be set
        #to an iterator over all changes that are to happen to the 
        #changing object during the model (where changes in the iterator
        #are tuples of (timestep, change_function))
        self.changes = None
        #create self.next_change, which defaults to None but will be
        #set by self.set_changes and upkept by self.make_change
        self.next_change = None

    def _set_next_change(self):
        try:
            self.next_change = next(self.changes)
        except StopIteration:
            self.next_change = None

    def _make_change(self, t, additional_args, for_plotting=False,
        verbose=False):
        #if this is the right timestep for the next change
        if self.next_change is not None:
            if t == self.next_change[0] and for_plotting == False:
                if verbose:
                    print("\t**** Running the next change\t%s\n\n" % str(
                        self.next_change))
                #call the next_change function to make the change
                #NOTE: unpacking additonal_args dict into it, so that the
                #LandChanger can receive the current Landscape object, rather
                #than capturing the Landscape object within the closure, taking
                #up memory, or capturing a pointer to the object and risking
                #problems that I'm not foreseeing (a la https://docs.python.
                #org/3/faq/programming.html#why-do-lambdas-defined-in-a
                #-loop-with-different-values-all-return-the-same-result),
                #and so that SpeciesChanger can similarly get the current
                #Species object
                self.next_change[1](changer = self, **additional_args)
                    #NOTE: feeding self in so that some change functions
                    #that need to set Changer attributes (e.g. spp-changer
                    #demographic functions set the base_K raster at
                    #the beginning of a period of dem changes)
                    #have access to the Changer object
                #then load the next change after that
                self._set_next_change()
                #then recurse this function, so that multiple changes
                #will happen this timestep if need be
                self._make_change(t, additional_args)

    #add a change (in the form of a (timestep, change_fn) pair) at the
    #appropriate timestep
    def _add_change(self, change):
        changes = [self.next_change]
        more = True
        while more:
            try:
                changes.append(next(self.changes))
            except StopIteration:
                more = False
        before = [c[0] <= change[0] for c in self.changes]
        insert = before.index(False)
        changes = changes[:insert] + [change] + changes[insert:]
        self.next_change = changes[0]
        self.changes = iter(changes)


class _LandscapeChanger(_Changer):
    def __init__(self, land, land_change_params, mod):
        super(_LandscapeChanger, self).__init__(land_change_params)
        #set the type-label of the changer
        self.type = 'land'

        #create an empty dict to hold all the basic change info
        #(end_rast, timesteps) for each layer that will change (so that
        #other objects, such as movement_surfaces, can access this
        #information and use it to build their own change-series)
        self.change_info = {}

    #call self.set_changes() to set self.changes and self.next_change
        self._set_changes(land)

    #method to set the changes stipulated in land_change_params dict for
    #the land object
    def _set_changes(self, land):

        #get the time-series of layers for each lyr_num
        lyrs = {}
        surfs = {}
        for lyr_num in self.change_params.keys():
            #get the conglomerate lyr series, across all change events
            #parameterized for this Layer
            conglom_lyr_series = _make_conglom_lyr_series(land, lyr_num,
                self.change_params[lyr_num])
            #add the conglomerate lyr series to the lyrs object
            lyrs[lyr_num] = conglom_lyr_series
            #set the _LandscapeChanger.change_info attribute, so that it can be
            #grabbed later for building _move_surf,
            #_disp_surf, and other objects
            self.change_info[lyr_num] = {**self.change_params[lyr_num]}

        #get everything in chronological order (with lyr-change functions
        #coming before surf-change functions if they're in the same timestep)
        lyr_changes = []
        for lyr_num in lyrs.keys():
            lyr_changes.extend([(t, lyr_num,
                                   lyr) for t,lyr in lyrs[lyr_num]])
        lyr_changes = sorted(lyr_changes, key = lambda x: x[0])

        #make a list of change_fns for the changes in all_changes, and make
        #self.changes an iterator over that list
        change_fns = [(change[0],
            _get_lyr_change_fn(*change[1:])) for change in lyr_changes]
        change_fns = iter(change_fns)
        self.changes = change_fns
        #load the first change to happen into the self.next_change attribute
        self._set_next_change()


class _SpeciesChanger(_Changer):
    def __init__(self, spp, spp_change_params, land):
        super(_SpeciesChanger, self).__init__(spp_change_params)
        self.type = 'spp'

        #an attribute that is used by some dem-change fns, as a baseline
        #population size at the start of the demographic change event
        #against which later timesteps' changes are defined
        self.base_K = None

    #call self._set_changes() to set self.changes and self.next_change
        self._set_changes(spp, land)

    #method to set the base_K attribute to spp.K
    def _set_base_K(self, spp):
        self.base_K = spp.K

    #method to set the changes stipulated in params dict for the spp object
    def _set_changes(self, spp, land):
        #pull out the parts of the params
        try:
            dem_change_params = self.change_params.dem
        except Exception:
            dem_change_params = None
        try:
            parameter_change_params = self.change_params.life_hist
        except Exception:
            parameter_change_params = None
        #check if this spp has a _move_surf, and if there are land-changes
        #that affect its lyr
        move_surf_change_fns = []
        if (spp._move_surf is not None
            and land._changer is not None
            and spp._move_surf.lyr_num in land._changer.change_info.keys()):
            #if so, grab that lyr's land change params and create a
            #move_surf_series with them, to be 
            #added to this spp's change fns
            lc_move_surf_params = land._changer.change_info[
                                                    spp._move_surf.lyr_num]
            #create a time-series of movement surfaces 
            surf = spp._move_surf
            surf_series = _make_conductance_surface_series(land,
                start_lyr = land[surf.lyr_num], mixture = surf.mix,
                kappa = surf.kappa, approx_len = surf.approx_len,
                change_params_one_lyr = lc_move_surf_params)
            #and create change fns from them
            move_surf_change_fns.extend(
                _get_conductance_surface_change_fns(surf_series))

        #check if this spp has a _disp_surf, and if there are land-changes
        #that affect its lyr
        disp_surf_change_fns = []
        if (spp._disp_surf is not None
            and land._changer is not None
            and spp._disp_surf.lyr_num in land._changer.change_info.keys()):
            #if so, grab that lyr's land change params and create a
            #disp_surf_series with them, to be 
            #added to this spp's change fns
            lc_disp_surf_params = land._changer.change_info[
                                                    spp._disp_surf.lyr_num]
            #create a time-series of movement surfaces 
            surf = spp._disp_surf
            surf_series = _make_conductance_surface_series(land,
                start_lyr = land[surf.lyr_num], mixture = surf.mix,
                kappa = surf.kappa, approx_len = surf.approx_len,
                change_params_one_lyr = lc_disp_surf_params)
            #and create change fns from them
            disp_surf_change_fns.extend(
                _get_conductance_surface_change_fns(surf_series))


        #set the demographic changes, if applicable
        dem_change_fns = []
        if dem_change_params is not None:
            for event, event_params in dem_change_params.items():
                if True in [v is not None for v in event_params.values()]:
                    dem_change_fns.extend(_get_dem_change_fns(spp,
                                                             **event_params))

        #set the other changes, if applicable
        parameter_change_fns = []
        if parameter_change_params is not None:
            for parameter, parameter_params in parameter_change_params.items():
                if True in [v is not None for v in parameter_params.values()]:
                    parameter_change_fns.extend(
                        _get_parameter_change_fns(parameter,
                        **parameter_params))

        #put all the changes in chronological order
        if (len(dem_change_fns) + len( parameter_change_fns) + len(
            move_surf_change_fns) + len(disp_surf_change_fns) > 0):
            change_fns = dem_change_fns + parameter_change_fns
            change_fns = change_fns + move_surf_change_fns
            change_fns = change_fns + disp_surf_change_fns
            change_fns = sorted(change_fns, key = lambda x: x[0])

            #make self.changes an iterator over that list
            change_fns = iter(change_fns)
            self.changes = change_fns
            #load the first change to happen into the self.next_change
            #attribute
            self._set_next_change()

    #a method to visualize the species changes that will occur
    def _plot_dem_changes(self, spp):
        if self.next_change is None:
            print(("No demographic changes remaining for this species. "
                "They were likely already run."))
        else:
            cop_spp = deepcopy(spp)
            cop_self = deepcopy(self)
            cop_changes = deepcopy(cop_self.changes)
            step_list = [cop_self.next_change[0]]
            more = True
            while more:
                try:
                    next_change = next(cop_changes)
                    step_list.append(next_change[0])
                except StopIteration:
                    more = False
            end = int(1.1*max(step_list))
            #set spp.K to 1, for viz purposes
            spp.K = 1
            #and set cop_self.base_K 
            cop_self._set_base_K(spp)
            Ks = []
            for t in range(end):
                cop_self._make_change(t, for_plotting = True,
                                      additional_args={'spp': spp})
                Ks.append(spp.K)
            plt.plot(range(end), Ks)
            #set spp back to its original value
            spp = cop_spp


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

    #########################
    # for _LandscapeChanger #
    #########################

# function that takes a starting rast, its layer, an ending rast, a few
# timestep arguments, and returns a linearly interpolated stack of rasters
def _make_lyr_series(lyr, change_rast, start_t, end_t, n_steps,
                     coord_prec=0):
    start_rast = lyr.rast
    dim = lyr.dim
    scale_min = lyr._scale_min
    scale_max = lyr._scale_max

    #get (rounded) evenly spaced timesteps at which to implement the changes
    timesteps = np.int64(np.round(np.linspace(start_t, end_t, n_steps)))

    #if change_rast is a file or an array
    if ((isinstance(change_rast, str) and os.path.isfile(change_rast))
        or isinstance(change_rast, np.ndarray)):
        #if it's an array, keep it, and set the prj, res, dim, ulc values;
        #if it's a file, read it in and get the dim, res, ulc, and prj values
        if isinstance(change_rast, np.ndarray):
            #NOTE: this resets dim from the argument that was fed in
            dim = change_rast.shape
            # NOTE: setting ulc, res, and prj to the same values as the layer
            # provided to the layer argument, so that the None values don't
            # throw an error in line 515's assertion below, but in reality
            # there might be a better, 'safer' way to handle this
            ulc = lyr.ulc
            res = lyr.res
            prj = lyr.prj
        elif os.path.isfile(change_rast):
            #NOTE: this resets dim from the argument that was fed in
            change_rast, dim, res, ulc, prj = io._read_raster(change_rast,
                                                              coord_prec,
                                                              dim=dim)
            # scale the raster correctly
            change_rast, scale_min_out, scale_max_out = spt._scale_raster(
                                            change_rast, scale_min, scale_max)
            assert (scale_min == scale_min_out
                    and scale_max == scale_max_out), ("The scale_min and "
                                                      "scale_max values "
                                                      "returned from scaling "
                                                      "the input change "
                                                      "raster don't match "
                                                      "those of Layer, to "
                                                      "which the raster "
                                                      "corresponds.") % lyr.idx
        #flatten the start and end rasters
        start = start_rast.flatten()
        end = change_rast.flatten()
        #get a column of values for each grid cell on the landscape, linearly
        #spaced between that cell's start and end values
        #NOTE: linspace(..., ..., n+1)[1:] gives us the changed rasters
        #for n steps, leaving off the starting lyr value because that's
        #already the existing lyr so we don't want that added into our changes
        rast_series = np.vstack([np.linspace(start[i],
            end[i], n_steps+1)[1:] for i in range(len(start))])
        #then reshape each timeslice into the dimensions of start_rast
        rast_series = [rast_series[:,i].reshape(
            start_rast.shape) for i in range(rast_series.shape[1])]

    #elif it's a directory
    elif isinstance(change_rast, str) and os.path.isdir(change_rast):
        files = os.listdir(change_rast)
        assert len(files) == n_steps, ("The number of files in the directory "
            "provided for the 'change_rast' parameter is not equal to the "
            "number provided for the 'n_steps' parameter.")
        #check that each file starts with an integer
        for f in files:
            assert os.path.splitext(f.split('_')[0])[0].isnumeric(), ("In the "
                "directory provided for the "
                "'change_rast' parameter, the file %s does not start "
                "with an integer followed by an underscore. This is "
                "necessary, to indicate the timestep during the Landscape "
                "change event at which the file's raster should "
                "be used.") % f
        #ensure that the files are ordered by the beginning integer
        #NOTE: should already be sorted by natural sorting, but just in case
        steps_and_files = {int(
            os.path.splitext(f.split('_')[0])[0]):f for f in files}
        files = [steps_and_files[i] for i in sorted(steps_and_files.keys())]
        #read in the whole series, check all prj, dim, ulc, and res values
        #are equal, then set the prj, dim, ulc, and res values, and
        #set rasters equal to rast_series
        all_rasts = [io._read_raster(os.path.join(
            change_rast, f), coord_prec, dim) for f in files]
        dim_cts = C([i[1] for i in all_rasts])
        assert len(dim_cts) == 1, ("The dimensions of "
            "all files in the directory provided for the 'change_rast' "
            "parameter are not equal. Most files have dimensions %s, but the "
            "following files differ: %s.") % (
            str([k for k,v in dim_cts.items() if v == max(dim_cts.values())]),
            str([files[k] for k,v in dim_cts.items() if v != max(
            dim_cts.values())]))
        #NOTE: this resets dim from the argument that was fed in
        dim = all_rasts[0][1]
        res_cts = C([tuple([*i[2]]) for i in all_rasts])
        assert len(res_cts) == 1, ("The spatial resolutions of "
            "all files in the directory provided for the 'change_rast' "
            "parameter are not equal. Most files have resolution %s, but the "
            "following files differ: %s.") % (
            str([k for k,v in res_cts.items() if v == max(res_cts.values())]),
            str([files[k] for k,v in res_cts.items() if v != max(
            res_cts.values())]))
        res = all_rasts[0][2]
        ulc_cts = C([tuple([*i[3]]) for i in all_rasts])
        assert len(ulc_cts) == 1, ("The upper left corners of "
            "all files in the directory provided for the 'change_rast' "
            "parameter are not equal. Most files have upper left corner %s, "
            "but the following files differ: %s.") % (
            str([k for k,v in ulc_cts.items() if v == max(ulc_cts.values())]),
            str([files[k] for k,v in ulc_cts.items() if v != max(
            ulc_cts.values())]))
        ulc = all_rasts[0][3]
        prj_cts = C([i[4] for i in all_rasts])
        assert len(prj_cts) == 1, ("The projections of "
            "all files in the directory provided for the 'change_rast' "
            "parameter are not equal. Most files have projection %s, but the "
            "following files differ: %s.") % (
            str([k for k,v in prj_cts.items() if v == max(prj_cts.values())]),
            str([files[k] for k,v in prj_cts.items() if v != max(
            prj_cts.values())]))
        prj = all_rasts[0][4]
        rast_series = [i[0] for i in all_rasts]
        # rescale the rasters
        rasts_and_scale_vals = [spt._scale_raster(
                        rast, scale_min, scale_max) for rast in rast_series]
        assert len(set([i[1] for i in rasts_and_scale_vals])) == 1, (
            "More than one unique scale_min value was returned from rescaling "
            "the directory of rasters provided for Layer %i.") % lyr.idx
        assert len(set([i[2] for i in rasts_and_scale_vals])) == 1, (
            "More than one unique scale_max value was returned from rescaling "
            "the directory of rasters provided for Layer %i.") % lyr.idx
        scale_min_out = rasts_and_scale_vals[0][1]
        scale_max_out = rasts_and_scale_vals[0][2]
        assert (scale_min == scale_min_out
                and scale_max == scale_max_out), ("The scale_min and "
                                                  "scale_max values "
                                                  "returned from scaling "
                                                  "the input change "
                                                  "rasters don't match "
                                                  "those of Layer %i, to "
                                                  "which the raster "
                                                  "corresponds.") % lyr.idx
        rast_series = [i[0] for i in rasts_and_scale_vals]
        #try to get the timesteps from the filenames
        try:
            timesteps = [int(os.path.splitext(
                i.split('_')[0])[0]) for i in files]
        except Exception as e:
            raise ValueError("Unable to extract timesteps from the beginning "
                "of all filenames in the directory provided to the "
                "'change_rast' parameter. Each filename must begin with an "
                "integer indicating the timestep at which that raster should "
                "be switched out for the Layer's previous value, followed "
                "by an underscore (e.g. '50_mat_2011.tif' for a raster that "
                "should be used at timestep 50).")
        #check that there are no repeat timesteps
        timestep_cts = C(timesteps)
        assert len(set(timestep_cts.values())) == 1, ("Not all timesteps "
            "indicated by the beginnings of the files in the directory "
            "provided for the 'change_rast' parameter are unique. "
            "Duplicated timesteps include: %s") % (
            str([k for k, v in timestep_cts if v > 1]))
        #check that first and last timesteps of the files equal start_t & end_t
        assert timesteps[0] == start_t, ("The timesteps indicated by the "
            "beginnings of the filenames must start at the timestep indicated "
            "by the 'start_t' parameter. The 'start_t' parameter was provided "
            "the value %i, but the first timestep found in the directory was "
            "%i (from the file named '%s').") % (start_t, timesteps[0],
            files[0])
        assert timesteps[::-1][0] == end_t, ("The timesteps indicated by the "
            "beginnings of the filenames must end at the timestep indicated "
            "by the 'end_t' parameter. The 'end_t' parameter was provided "
            "the value %i, but the last timestep found in the directory was "
            "%i (from the file named '%s').") % (end_t, timesteps[::-1][0],
            files[::-1][0])

    #otherwise raise a ValueError
    else:
        raise ValueError(("The value provided for the 'change_rast' parameter "
            "must be either a numpy.ndarray or a path to a valid file "
            "(indicating the endpoint-raster of the change event), or a "
            "path to a directory of valid files (one for each timestep in "
            "the change event)."))

    #check that all the lengths match up
    assert len(rast_series) == n_steps, ("The length of the rast_series "
        "variable is not equal to the n_steps variable.")
    assert len(rast_series) == len(timesteps), ("The number of changing "
        "rasters created is not the same as the number of timesteps to be "
        "assigned to them")

    #zip the timesteps and the rasters together and return them as a list
    rast_series = list(zip(timesteps, rast_series))

    return(rast_series, dim, res, ulc, prj)


#function for getting the congolmerate lyr_series for all change
#events parameterized for a Layer
def _make_conglom_lyr_series(land, lyr_num, change_params_one_lyr):
    #check that the timeframes of the change events for this layer
    #don't overlap 
    all_timesteps = [[*range(v.start_t,
        v.end_t)] for v in change_params_one_lyr.values()]
    all_timesteps = [t for i in all_timesteps for t in i]
    assert (len(set(all_timesteps)) == len(all_timesteps)), ("Some "
        "of the change events for Layer number %i overlap "
        "in time.") % lyr_num
    #get the congolmerated lyr_series for all change events for this
    #layer
    all_lyr_series = [_make_lyr_series(land[lyr_num],
                                       coord_prec=land[lyr_num].coord_prec,
                                **v) for v in change_params_one_lyr.values()]
    #check that dim, res, ulc, and prj values are the same for each
    #change event
    assert len(set([i[1] for i in all_lyr_series]))==1, (""
        "Not all change events for Layer number %i produce layer "
        "series with the same dimensions.  Resulting dimensions "
        "include: (%s).") % (lyr_num, ', '.join([*set(
        [i[1] for i in all_lyr_series])]))
    assert len(set([tuple([*i[2]]) for i in all_lyr_series]))==1, (""
        "Not all change events for Layer number %i produce layer "
        "series with the same resolution.  Resulting resolutions "
        "include: (%s).") % (lyr_num, ', '.join([*set(
        [i[2] for i in all_lyr_series])]))
    assert len(set([tuple([*i[3]]) for i in all_lyr_series]))==1, (""
        "Not all change events for Layer number %i produce layer "
        "series with the same upper-left corners.  Resulting "
        "upper-left corners include: (%s).") % (lyr_num, ', '.join([*set(
        [i[3] for i in all_lyr_series])]))
    assert len(set([i[4] for i in all_lyr_series]))==1, (""
        "Not all change events for Layer number %i produce layer "
        "series with the same projection.  Resulting projections "
        "include: (%s).") % (lyr_num, ', '.join([*set(
        [i[4] for i in all_lyr_series])]))
    #take dim, res, ulc, and prj as single values
    dim_res_ulc_prj_list = [[i[n] for i in all_lyr_series] for n in range(1,5)]
    # NOTE: next line just coerces ndarrays to tuples, so that they can be
    # hashed into a set in the following line
    dim_res_ulc_prj_list = [[tuple([*i]) for i in item] if isinstance(
        item[0], np.ndarray) else item for item in dim_res_ulc_prj_list]
    dim, res, ulc, prj = [[*set(item)][0] for item in dim_res_ulc_prj_list]
    #check that dim, res, ulc, and prj match that of the land
    assert dim == land.dim, ('Dimensionality of lyr_series fed '
        'into _LandscapeChanger for Layer number %i does not '
        'match that of the land to be changed.') % lyr_num
    assert res == land.res or res is None, ('Resolution of '
        'lyr_series fed into _LandscapeChanger for Layer number %i '
        'does not match that of the land to be changed.') % lyr_num
    assert ulc == land.ulc or ulc is None, ('Upper-left corner of '
        'lyr_series fed into _LandscapeChanger for Layer number %i '
        'does not match that of the land to be changed.') % lyr_num
    assert prj == land.prj or prj is None, ('Projection of '
        'lyr_series fed into _LandscapeChanger for Layer number %i '
        'does not match that of the land to be changed.') % lyr_num
    #combine all lyr_series into one single lyr series
    conglom_lyr_series = [lyr_step for series in [
        i[0] for i in all_lyr_series] for lyr_step in series]
    return conglom_lyr_series


#function that takes a Landscape, a lyr_num and a dictionary of
#{t_change:new_lyr} and creates an _LandscapeChanger object that will change
#out Landscape[lyr_num] for new_lyr at each requisite t_change timestep
def _make_custom_lyr_series(lyr_dict):
    pass


def _get_lyr_change_fn(lyr_num, new_lyr_rast):
    def fn(changer, land, lyr_num = lyr_num, new_lyr_rast = new_lyr_rast):
        land._set_raster(lyr_num, new_lyr_rast)
    return(fn)


    #######################
    # for _SpeciesChanger #
    #######################
def _make_conductance_surface_series(land, start_lyr, mixture, kappa,
        approx_len, change_params_one_lyr):
    #get the time-series of lyrs across the land-change event (can
    #be at a reduced temporal resolution determined by t_res_reduct_factor,
    #because probably unnecessary to change the move_surf or disp_surf every
    #time the land changes, and they could be a fairly large objects to
    #hold in memory anyhow; but t_res_reduct_factor defaults to 1)
    conglom_lyr_series= _make_conglom_lyr_series(land, start_lyr.idx,
        change_params_one_lyr)
    #then get the series of _ConductanceSurface objects
    surf_series = []
    dummy_lyr = deepcopy(start_lyr)
    for t, rast in conglom_lyr_series:
        #change the lyr's raster
        dummy_lyr.rast = rast
        #then create a Movement_Surface using the copied Landscape
        surf_series.append((t, spt._ConductanceSurface(dummy_lyr, mixture,
            approx_len = approx_len, vm_distr_kappa = kappa)))
    return(surf_series)


def _get_conductance_surface_change_fns(surf_series):
    timesteps = []
    fns = []
    for t, new_surf in surf_series:
        def fn(changer, spp, new_surf = new_surf):
            spp._move_surf = new_surf
        timesteps.append(t)
        fns.append(fn)
    change_fns = zip(timesteps, fns)
    return(change_fns)

        ###############
        # demographic #
        ###############

def _get_dem_change_fns(spp, kind, start_t=None, end_t=None, rate=None,
    interval=None, n_cycles=None, size_range=None,
    distr='uniform', min_size=None, max_size=None, timesteps=None,
                                    sizes=None, increase_first=True):
    if kind == 'monotonic':
        fns = _get_monotonic_dem_change_fns(rate = rate,
            start_t = start_t, end_t = end_t)
    elif kind == 'stochastic':
        fns = _get_stochastic_dem_change_fns(start_K = spp.K,
            start_t = start_t, end_t = end_t, interval = interval,
            size_range = size_range, distr = distr)
    elif kind == 'cyclical':
        fns = _get_cyclical_dem_change_fns(start_t = start_t,
            end_t = end_t, n_cycles = n_cycles, size_range = size_range,
            min_size = min_size, max_size = max_size,
                                        increase_first = increase_first)
    elif kind == 'custom':
        fns = _get_custom_dem_change_fns(timesteps = timesteps, sizes = sizes)
    return(fns)


def _make_dem_change_fns(sizes, timesteps, K_mode='base'):
    fns = []
    if K_mode == 'current':
        for size in sizes:
            def fn(changer, spp, size = size):
                spp.K*=size
            fns.append(fn)
    elif K_mode == 'base':
        t0 = timesteps[0]
        for size in sizes:
            def fn(changer, spp, size=size, t0 = t0):
                if spp.t == t0:
                    changer._set_base_K(spp)
                spp.K = changer.base_K*size
            fns.append(fn)
    change_fns = list(zip(timesteps, fns))
    return(change_fns)


#will generate exponential change in the population size by iteratively
#multiplying a carrying-capacity (K) raster
def _get_monotonic_dem_change_fns(rate, start_t, end_t):
    #get the timesteps for the demogprahic changes
    #NOTE: setting start_t and end_t to the same value will create a
    #single-timestep change (e.g. a sudden bottleneck or rapid expansion)
    timesteps = range(start_t, end_t+1)
    sizes = [rate]*len(timesteps)
    change_fns = _make_dem_change_fns(sizes, timesteps, K_mode='current')
    return(change_fns)


#NOTE: should I provide an option for linear monotonic change (because I can
#set pc.base_K and then multiply it by rate*t at each t of a period of change)?


#will create stochastic changes around a baseline carrying-capacity (K) raster
def _get_stochastic_dem_change_fns(start_K, size_range, start_t, end_t,
    interval, distr = 'uniform'):
    if interval is None:
        interval = 1
    timesteps = range(start_t, end_t + 1, interval)
    if distr == 'uniform':
        sizes = r.uniform(*size_range, len(timesteps))
    elif distr == 'normal':
        mean = np.mean(size_range)
        sd = (size_range[1] - size_range[0])/6
        sizes = r.normal(loc = mean, scale = sd, size = len(timesteps))
    else:
        raise ValueError(("Argument 'distr' must be a value among "
                                            "['uniform', 'normal']"))
    #make size return to starting size
    sizes[-1] = 1
    #make all the change functions
    change_fns = _make_dem_change_fns(sizes, timesteps, K_mode = 'base')
    return(change_fns)


def _get_cyclical_dem_change_fns(start_t, end_t, n_cycles,
    size_range=None, min_size=None, max_size=None, increase_first=True):
    #detemine the min and max sizes for the cycles, based on input arguments
    if size_range is not None and min_size is None and max_size is None:
        min_size, max_size = size_range
    elif size_range is None and min_size is not None and max_size is not None:
        pass
    #and throw an informative error if both size_range and min_size&max_size
    #arguments are provided (or neither)
    else:
        raise ValueError(('Must either provide size_range (as a tuple of '
            'minimum and maximum sizes), or provide min_size and max_size '
                                                'separately, but not both.'))
    #check that there are enough timesteps to complete the cycles
    assert n_cycles <= (end_t - start_t)/2, ('The number of cycles requested must '
        'be no more than half the number of time steps over which the cycling '
                                                        'should take place.')
    #create a base cycle (one sine cycle)
    base = np.sin(np.linspace(0, 2*np.pi,1000))
    #flip it if cycles should decrease before increasing
    if not increase_first:
        base = base[::-1]
    #scale the base cycle to the min_size
    scaled_base = [1 + n*(max_size-1) if n>=0 else n for n in base]
    #and scale it to the max_size
    scaled_base = [1 + n*(1-min_size) if n<0 else n for n in scaled_base]
    scaled_base = np.array(scaled_base)
    #then get the timesteps at which each cycle should complete
    #NOTE: subtract 1 to imitate Python range stop-exclusivity
    cycle_timesteps = np.int32(np.linspace(start_t, end_t, n_cycles+1))
    #get from that the length in timesteps of each cycle
    cycle_lengths = np.diff(cycle_timesteps)
    #get from that the pop sizes at each timestep (which will approximate
    #the scaled sine curve from scaled_base
    sizes = np.hstack([scaled_base[np.int32(np.linspace(1, len(scaled_base)-1,
                                            l))] for l in cycle_lengths] + [1])
    #then get all timesteps across the cycles, to match up to the pop sizes
    timesteps = range(cycle_timesteps[0], cycle_timesteps[-1]+1)
    #get all the change functions for those timesteps
    change_fns = _make_dem_change_fns(sizes, timesteps, K_mode = 'base')
    return(change_fns)


def _get_custom_dem_change_fns(timesteps, sizes):
    assert len(timesteps) == len(sizes), ('For custom demographic changes, '
                    'timesteps and sizes must be iterables of equal length.')
    change_fns = _make_dem_change_fns(sizes, timesteps, K_mode = 'base')
    return(change_fns)


        #################
        #   life_hist   #
        #################

def _make_parameter_change_fns(parameter, timesteps, vals):
    fns = []
    for val in vals:
        def fn(changer, spp, parameter = parameter, val = val):
            setattr(spp, parameter, val)
        fns.append(fn)
    change_fns = list(zip(timesteps, fns))
    return(change_fns)


def _get_parameter_change_fns(parameter, timesteps, vals):
    assert len(timesteps) == len(vals), ("For custom changes of the '%s' "
        "parameter, timesteps and vals must be iterables of equal "
        "length.") % parameter
    change_fns = _make_parameter_change_fns(parameter, timesteps, vals)
    return(change_fns)

