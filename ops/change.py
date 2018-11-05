#!/usr/bin/python
# change.py

'''
##########################################

Module name:          sims.change


Module contains:
                      - classes and functions to facilitate environmental,
                        demographic, and parameter change


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           07-06-18
Documentation:        URL


##########################################
'''
#geonomics imports
from utils import spatial as spt

#other imports
import numpy as np
import numpy.random as r
import matplotlib.pyplot as plt
from operator import attrgetter as ag
from collections import OrderedDict as OD
from copy import deepcopy
import itertools as it
from copy import deepcopy


#TODO:

    #1 I created a quick, crude way of stipulating other-parameter pop changes,
    #but it would probably be nicer to just further generalize the
    #dem-change functions to allow for linear, stochastic, cyclic, or custom
    #changes of other param values?

    #2 I think this can be standardized (e.g. _LandscapeChanger
    #get_<   >_change_fns functions I believe return just the
    #functions, whereas _PopulationChanger
    #ones return list(zip(t, change_fn)) objects) and generalized still
    #further. Worth doing now, while the code is fresh in my head...


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

#base class for _LandscapeChanger and _PopulationChanger
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

    def _make_change(self, t, for_plotting=False):
        #if this is the right timestep for the next change
        if self.next_change is not None:
            if t == self.next_change[0] and for_plotting == False:
                print("\t**** Running the next change\t%s\n\n" % str(
                                                    self.next_change))
                #call the next_change function to make the change
                self.next_change[1](self)
                    #NOTE: feeding self in so that some change functions
                    #that need to set Changer attributes (e.g. pop-changer
                    #dem functions set the #base_K raster at the beginning
                    #of a period of dem changes) have access to the
                    #Changer object
                #then load the next change after that
                self._set_next_change()
                #then recurse this function, so that multiple changes
                #will happen this timestep if need be
                self._make_change(t)

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
    def __init__(self, land, land_change_params):
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

        #get the linearly spaced time-series of layer for each lyr_num
        lyrs = {}
        surfs = {}
        for lyr_num in self.change_params.keys():
            lyr_series, dim, res, ulc, prj = _make_linear_lyr_series(
                        land[lyr_num].rast, **self.change_params[lyr_num])
            assert dim == land.dim, ('Dimensionality of lyr_series fed into '
                '_LandscapeChanger does not match that of the land '
                'to be changed.')
            assert res == land.res or res is None, ('Resolution of '
                'lyr_series fed into _LandscapeChanger does not match '
                'that of the land to be changed.')
            assert ulc == land.ulc or ulc is None, ('Upper-left corner of '
                'lyr_series fed into _LandscapeChanger does not match that of '
                'the land to be changed.')
            assert prj == land.prj or prj is None, ('Projection of '
            'lyr_series fed into _LandscapeChanger does not match that of the '
            'land to be changed.')
            lyrs[lyr_num] = lyr_series

            #set the _LandscapeChanger.change_info attribute, so that it can be
            #grabbed later for building _move_surf,
            #_disp_surf, and other objects
            self.change_info[lyr_num] = {**self.change_params[lyr_num]}
            self.change_info[lyr_num]['end_rast'] = lyr_series[-1][1]

        #get everything in chronological order (with lyr-change functions
        #coming before surf-change functions if they're in the same timestep)
        lyr_changes = []
        for lyr_num in lyrs.keys():
            lyr_changes.extend([(t, lyr_num,
                                   lyr) for t,lyr in lyrs[lyr_num]])
        lyr_changes = sorted(lyr_changes, key = lambda x: x[0])

        #make a list of change_fns for the changes in all_changes, and make
        #self.changes an iterator over that list
        change_fns = [(change[0], _get_lyr_change_fn(land,
                                    *change[1:])) for change in lyr_changes]
        change_fns = iter(change_fns)
        self.changes = change_fns
        #load the first change to happen into the self.next_change attribute
        self._set_next_change()


class _PopulationChanger(_Changer):
    def __init__(self, pop, pop_change_params, land=None):
        super(_PopulationChanger, self).__init__(pop_change_params)
        self.type = 'pop'

        #an attribute that is used by some dem-change fns, as a baseline
        #population size at the start of the demographic change event
        #against which later timesteps' changes are defined
        self.base_K = None

    #call self._set_changes() to set self.changes and self.next_change
        self._set_changes(pop)

    #method to set the base_K attribute to pop.K
    def _set_base_K(self, pop):
        self.base_K = pop.K

    #method to set the changes stipulated in params dict for the pop object
    def _set_changes(self, pop, land=None):
        #pull out the parts of the params
        try:
            dem_change_params = self.change_params.dem
        except Exception:
            dem_change_params = None
        try:
            parameter_change_params = self.change_params.life_hist
        except Exception:
            parameter_change_params = None
        #check if this pop has a _move_surf, and if there are land-changes
        #that affect its lyr
        move_surf_change_fns = []
        if (pop._move_surf is not None
            and land is not None
            and pop._move_surf.lyr_num in land._changer.change_info.keys()):
            #if so, grab that lyr's land change params and create a
            #move_surf_series with them, to be 
            #added to this pop's change fns
            lc_move_surf_params = land._changer.change_info[
                                                    pop._move_surf.lyr_num]
            #create a time-series of movement surfaces 
            surf = pop._move_surf
            surf_series = _make_conductance_surface_series(
                start_lyr = land[surf.lyr_num], mixture = surf.mix,
                kappa = surf.kappa, approx_len = surf.approx_len,
                **lc_move_surf_params)
            #and create change fns from them
            move_surf_change_fns.extend(_get_conductance_surface_change_fns(
                                                        pop, surf_series))

        #check if this pop has a _disp_surf, and if there are land-changes
        #that affect its lyr
        disp_surf_change_fns = []
        if (pop._disp_surf is not None
            and land is not None
            and pop._disp_surf.lyr_num in land._changer.change_info.keys()):
            #if so, grab that lyr's land change params and create a
            #disp_surf_series with them, to be 
            #added to this pop's change fns
            lc_disp_surf_params = land._changer.change_info[
                                                    pop._disp_surf.lyr_num]
            #create a time-series of movement surfaces 
            surf = pop._disp_surf
            surf_series = _make_conductance_surface_series(
                start_lyr = land[surf.lyr_num], mixture = surf.mix,
                kappa = surf.kappa, approx_len = surf.approx_len,
                **lc_disp_surf_params)
            #and create change fns from them
            disp_surf_change_fns.extend(_get_conductance_surface_change_fns(
                                                        pop, surf_series))


        #set the demographic changes, if applicable
        dem_change_fns = []
        if dem_change_params is not None:
            for event, event_params in dem_change_params.items():
                if True in [v is not None for v in event_params.values()]:
                    dem_change_fns.extend(_get_dem_change_fns(pop,
                                                             **event_params))

        #set the other changes, if applicable
        parameter_change_fns = []
        if parameter_change_params is not None:
            for parameter, parameter_params in parameter_change_params.items():
                if True in [v is not None for v in parameter_params.values()]:
                    parameter_change_fns.extend(_get_parameter_change_fns(pop,
                                        parameter, **parameter_change_params))

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

    #a method to visualize the population changes that will occur
    def _plot_dem_changes(self, pop):
        if self.next_change is None:
            print(("No demographic changes remaining for this population. "
                "They were likely already run."))
        else:
            cop_pop = deepcopy(pop)
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
            #set pop.K to 1, for viz purposes
            pop.K = 1
            #and set cop_self.base_K 
            cop_self._set_base_K(pop)
            Ks = []
            for t in range(end):
                cop_self._make_change(t, for_plotting = True)
                Ks.append(pop.K)
            plt.plot(range(end), Ks)
            #set pop back to its original value
            pop = cop_pop


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

    #########################
    # for _LandscapeChanger #
    #########################

#function that takes a starting lyr, an ending lyr, a number of
#timesteps, and a Model object, and returns a linearly interpolated
#stack of rasters
def _make_linear_lyr_series(start_rast, end_rast, start_t, end_t, n_steps):
    assert start_rast.shape == end_rast.shape, ('The starting raster and '
        'ending raster for the land-change event are not of the same '
        'dimensions: START: %s,  END %s') % (str(start_rast.shape),
                                                    str(end_rast.shape))

    if type(end_rast) is str:
        end_rast, dim, res, ulc, prj = io._read_raster(end_rast)

    elif type(end_rast) is np.ndarray:
        dim = end_rast.shape
        ulc = None
        res = None
        prj = None

    #get (rounded) evenly spaced timesteps at which to implement the changes
    timesteps = np.int64(np.round(np.linspace(start_t, end_t, n_steps)))
    #flatten the starting and ending rasters
    start = start_rast.flatten()
    end = end_rast.flatten()
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
    #check that all the lenghts match up
    assert len(rast_series) == n_steps, ("The length of the rast_series "
        "variable is not equal to the n_steps variable.")
    assert len(rast_series) == len(timesteps), ("The number of changing "
        "rasters created is not the same as the number of timesteps to be "
        "assigned to them")
    #zip the timesteps and the rasters together and return them as a list
    rast_series = list(zip(timesteps, rast_series))
    return(rast_series, dim, res, ulc, prj)


#function that takes a Landscape, a lyr_num and a dictionary of
#{t_change:new_lyr} and creates an _LandscapeChanger object that will change
#out Landscape[lyr_num] for new_lyr at each requisite t_change timestep
def _make_custom_lyr_series(lyr_dict):
    pass


def _get_lyr_change_fn(land, lyr_num, new_lyr):
    def fn(lc, land = land, lyr_num = lyr_num, new_lyr = new_lyr):
        land._set_raster(lyr_num, new_lyr)
    return(fn)


    ##########################
    # for _PopulationChanger #
    ##########################

def _make_conductance_surface_series(start_lyr, mixture, kappa,
        approx_len, end_rast, start_t, end_t, n_steps,
        t_res_reduct_factor=1):
    #get the time-series of lyrs across the land-change event (can
    #be at a reduced temporal resolution determined by t_res_reduct_factor,
    #because probably unnecessary to change the move_surf or disp_surf every
    #time the land changes, and they could be a fairly large objects to
    #hold in memory anyhow; but t_res_reduct_factor defaults to 1)
    lyr_series, dim, res, ulc, prj = _make_linear_lyr_series(
        start_lyr.rast, end_rast, start_t, end_t, int(
                                    n_steps*t_res_reduct_factor))
    #then get the series of _ConductanceSurface objects
    surf_series = []
    dummy_lyr = deepcopy(start_lyr)
    for t, rast in lyr_series:
        #change the lyr's raster
        dummy_lyr.rast = rast
        #then create a Movement_Surface using the copied Landscape
        surf_series.append((t, spt._ConductanceSurface(dummy_lyr, mixture,
            approx_len = approx_len, vm_distr_kappa = kappa)))
    return(surf_series)


def _get_conductance_surface_change_fns(pop, surf_series):
    timesteps = []
    fns = []
    for t, new_surf in surf_series:
        def fn(pc, pop=pop, new_surf = new_surf):
            pop._move_surf = new_surf
        timesteps.append(t)
        fns.append(fn)
    change_fns = zip(timesteps, fns)
    return(change_fns)

        ###############
        # demographic #
        ###############

def _get_dem_change_fns(pop, kind, start=None, end=None, rate=None,
    interval=None, n_cycles=None, size_range=None,
    distr='uniform', min_size=None, max_size=None, timesteps=None,
                                    sizes=None, increase_first=True):
    if kind == 'monotonic':
        fns = _get_monotonic_dem_change_fns(pop = pop, rate = rate,
                    start = start, end = end)
    elif kind == 'stochastic':
        fns = _get_stochastic_dem_change_fns(pop = pop, start = start,
            end = end, interval = interval, size_range = size_range,
                                                            distr = distr)
    elif kind == 'cyclical':
        fns = _get_cyclical_dem_change_fns(pop = pop, start = start,
            end = end, n_cycles = n_cycles, size_range = size_range,
            min_size = min_size, max_size = max_size,
                                        increase_first = increase_first)
    elif kind == 'custom':
        fns = _get_custom_dem_change_fns(pop = pop, timesteps = timesteps,
                                                            sizes = sizes)
    return(fns)


def _make_dem_change_fns(pop, sizes, timesteps, K_mode='base'):
    fns = []
    if K_mode == 'current':
        for size in sizes:
            def fn(pc, pop = pop, size = size):
                pop.K*=size
            fns.append(fn)
    elif K_mode == 'base':
        t0 = timesteps[0]
        for size in sizes:
            def fn(pc, pop=pop, size=size, t0 = t0):
                if pop.t == t0:
                    pc._set_base_K(pop)
                pop.K = pc.base_K*size
            fns.append(fn)
    change_fns = list(zip(timesteps, fns))
    return(change_fns)


#will generate exponential change in the population by iteratively
#multiplying a carrying-capacity (K) raster
def _get_monotonic_dem_change_fns(pop, rate, start, end):
    #get the timesteps for the demogprahic changes (subtract 1 from start
    #to set Python's 0th timestep as 1)
    #NOTE: setting start and end to the same value will create a
    #single-timestep change (e.g. a sudden bottleneck or rapid expansion)
    timesteps = range(start, end)
    sizes = [rate]*len(timesteps)
    change_fns = _make_dem_change_fns(pop, sizes, timesteps, K_mode='current')
    return(change_fns)


#NOTE: should I provide an option for linear monotonic change (because I can
#set pc.base_K and then multiply it by rate*t at each t of a period of change)?


#will create stochastic changes around a baseline carrying-capacity (K) raster
def _get_stochastic_dem_change_fns(pop, size_range, start, end, interval,
                                                        distr = 'uniform'):
    start_K = pop.K
    if interval is None:
        interval = 1
    timesteps = range(start, end, interval)
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
    change_fns = _make_dem_change_fns(pop, sizes, timesteps, K_mode = 'base')
    return(change_fns)


def _get_cyclical_dem_change_fns(pop, start, end, n_cycles, size_range=None,
                            min_size=None, max_size=None, increase_first=True):
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
    assert n_cycles <= (end - start)/2, ('The number of cycles requested must '
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
    cycle_timesteps = np.int32(np.linspace(start, end, n_cycles+1))
    #get from that the length in timesteps of each cycle
    cycle_lengths = np.diff(cycle_timesteps)
    #get from that the pop sizes at each timestep (which will approximate
    #the scaled sine curve from scaled_base
    sizes = np.hstack([scaled_base[np.int32(np.linspace(1, len(scaled_base)-1,
                                            l))] for l in cycle_lengths] + [1])
    #then get all timesteps across the cycles, to match up to the pop sizes
    timesteps = range(cycle_timesteps[0], cycle_timesteps[-1]+1)
    #get all the change functions for those timesteps
    change_fns = _make_dem_change_fns(pop, sizes, timesteps, K_mode = 'base')
    return(change_fns)


def _get_custom_dem_change_fns(pop, timesteps, sizes):
    assert len(timesteps) == len(sizes), ('For custom demographic changes, '
                    'timesteps and sizes must be iterables of equal length.')
    change_fns = _make_dem_change_fns(pop, sizes, timesteps, K_mode = 'base')
    return(change_fns)


        #################
        #   life_hist   #
        #################

def _make_parameter_change_fns(pop, parameter, timesteps, vals):
    fns = []
    for val in vals:
        def fn(pc, parameter, pop = pop, val = val):
            setattr(pop, parameter, val)
        fns.append(fn)
    change_fns = list(zip(timesteps, fns))
    return(change_fns)


def _get_parameter_change_fns(pop, parameter, timesteps, vals):
    assert len(timesteps) == len(vals), ("For custom changes of the '%s' "
     "paramter, timesteps and vals must be iterables of equal length.") % param
    change_fns = _make_parameter_change_fns(pop, param, timesteps, vals)
    return(change_fns)

