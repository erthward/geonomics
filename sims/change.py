#!/usr/bin/python
# change.py

'''
##########################################

Module name:          sims.change


Module contains:
                      - classes and functions to facilitate environmental, demographic, and parameter change


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           07-06-18
Documentation:        URL


##########################################
'''
#geonomics imports
from utils import spatial as spt
spt.Movement_Surface

#other imports
import numpy as np
from operator import attrgetter as ag
from collections import OrderedDict as OD
import itertools as it
import copy


#TODO:
    #3 rename movement_surf to move_surf package-wide?

    #4 then think about how to further generalize these structures, and create similar ones for Gen_Changer
    #and Pop_Changer

    #4.5 simplest way to express range of possible demographic scenarios in the params file, and then
    #interpret those and split those into different demographic-change functions here?

    #5 then write the overall Sim_Changer class to unite everything

    #6 try running and debugging what I've got so far, and reflect on if it's a decent wa to do this or if it
    #needs retooling



######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

#base class for Land_, Pop_, and Gen_Changer classes
class Changer:
    def __init__(self):
        #set the type-label of the changer
        self.type = None

        #create self.changes, which defaults to None but will be set to an iterator over all changes
        #that are to happen to the changing object during the model (where changes in the iterator are tuples 
        #of (timestep, change_function))
        self.changes = None
        #create self.next_change, which defaults to None but will be set by self.set_changes and upkept by self.make_change
        self.next_change = None

        #set self.changes and self.next_change
        #self.set_changes(land, params)
       
    def set_next_change(self):
        self.next_change = next(self.changes)

    def make_change(self, t):
        #if this is the right timestep for the next change
        if t == self.next_change[0]:
            print("\n\nRUNNING THE NEXT CHANGE:\n\t%s" % str(self.next_change))
            #call the next_change function to make the change
            self.next_change[1]()
            #then load the next change after that
            self.set_next_change()
            #then recurse this function, so that multiple changes will happen this timestep if need be
            self.make_change(t)


class Land_Changer(Changer):
    def __init__(self, land, params):
        super(Land_Changer, self).__init__()
        #set the type-label of the changer
        self.type = 'land'

    #call self.set_changes() to set self.changes and self.next_change
        self.set_changes(land, params)

    #method to set the changes stipulated in params dict for the land object
    def set_changes(self, land, params):
        #pull out the separate parts of the params file
        land_change_params = params['change']['land']

        #get the linearly spaced time-series of scapes for each scape_num
        scapes = {}
        surfs = {}
        for scape_num in land_change_params.keys():
            scape_series = make_linear_scape_series(land.scapes[scape_num].raster, **land_change_params[scape_num])
            scapes[scape_num] = scape_series
            if scape_num == land.movement_surf_scape_num:
                surf_series = make_movement_surf_series(land, scape_series)
                surfs[scape_num] = surf_series
        
        #check that timesteps are the same for the surf_series and scape_series items pertaining to the
        #movement-surface scape
        if land.movement_surf_scape_num in scapes.keys():
            assert [i[0] for i in scapes[land.movement_surf_scape_num]] == [i[0] for i in surfs[land.movement_surf_scape_num]], 'ERROR: scape_series and surf_series do not contain the same timesteps for the movement surface scape'
        
        #get everything in chronological order (with scape-change functions coming before surf-change
        #functions if they're in the same timestep)
        scape_changes = []
        for scape_num in scapes.keys():
            scape_changes.extend([(t, scape_num, scape) for t,scape in scapes[scape_num]])
            #if scape_num == land.movement_surf_scape_num:
                #all_changes.extend([(t, scape_num, surf) for t,surf in surfs[scape_num]])
        scape_changes = sorted(scape_changes, key = lambda x: x[0])
        all_changes = []
        for change in scape_changes:
            all_changes.append(change)
            if change[1] == land.movement_surf_scape_num:
                corresp_surf_change = [(t, change[1], surf) for t, surf in surfs[change[1]] if list(surfs.keys())[0] == change[1] and t == change[0]]
                assert len(corresp_surf_change) == 1, "ERROR: Found more than one surface-change function corresponding to a given scape-change function (i.e. for a given timestep and scape_num combo)"
                all_changes.append(corresp_surf_change[0])

        #make a list of change_fns for the changes in all_changes, and make self.changes an iterator over that list
        change_fns = [(change[0], get_land_change_fn(land, *change[1:])) for change in all_changes]
        change_fns = iter(change_fns)
        self.changes = change_fns
        #load the first change to happen into the self.next_change attribute
        self.set_next_change()


class Pop_Changer:
    def __init__(self, change_vals, params):
        self.type = 'pop'

    #NOTE: in the following pop_size can just always be expressed as a fraction of the starting pop (e.g 1.1
        #for a 1000-pop is 1100) or as absolute size, with a params-arg choosing which('frac', 'abs')
    #I think this should probably take dem arguments such as 'custom' (list of (t, pop_size) tuples), 
                                                            #'bottleneck' ('instant', or 'gradual')
                                                            #'decrease'   ( rate and period of exponential decrease, or target size)
                                                            #'increase'   ( rate and period of exponential increase, or target size)
                                                            #'stochastic' ( range of fluctuations)
                                                            #'cyclical'
    #NOTE: changes in the environment will also of course drive pop-size changes intrinsically... how to 
    #square this? does this really just boil down to Land_Changer changes on the hab/K raster?  
    
    #call self.set_changes() to set self.changes and self.next_change
        self.set_changes(land, params)

    #method to set the changes stipulated in params dict for the pop object
    def set_changes(self, pop, params):
        #pull out the relevant part of the params file
        pop_change_params = params['change']['pop']
        dem_change_params = pop_change_params['dem']
        other_change_params = pop_change_params['other']

        #set the demographic changes, if applicable
        if len(dem_change_params) > 0 and True in [v is not None for v in dem_change_params.values()]:
            dem_change_fns = get_dem_change_fns(**dem_change_params)

        #set the other changes, if applicable







class Gen_Changer:
    def __init__(self, change_vals, params):
        self.type = 'gen'

#an overall Sim_Changer class, which will contain and operate the Land_, Pop_, and Gen_Changer objects
#and will make each timestep's necessary changes in that order (land-level changes, then population, then genome)
class Sim_Changer:
    def __init__(self, changers, params):
        type_getter = ag('type')
        valid_changer_types = ['land', 'pop', 'gen']
        types = list(map(type_getter, changers))
        for t in types:
            setattr(self, t, [changer for changer in changers if changer.type == t])
        for t in set(valid_changer_types) - set(types):
            setattr(self, t, None)


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

    ####################
    # for Land_Changer #
    ####################
    
#function that takes a starting scape, an ending scape, a number of timesteps, and a Model object,
#and returns a linearly interpolated stack of scapes
def make_linear_scape_series(start_scape, end_scape, t_start, t_end, n):
    #get (rounded) evenly spaced timesteps at which to implement the changes
    timesteps = np.int64(np.round(np.linspace(t_start, t_end, n)))
    #flatten the starting and ending scapes
    start = start_scape.flatten()
    end = end_scape.flatten()
    #get a column of values for each grid cell on the landscape, linearly spaced between that cell's start and end values
    #NOTE: linspace(..., ..., n+1)[1:] gives us the changed scapes for n steps, leaving off the starting scape
    #value because that's already the existing scape so we don't want that added into our changes
    scape_series = np.vstack([np.linspace(start[i], end[i], n+1)[1:] for i in range(len(start))])
    #then reshape each timeslice into the dimensions of start_scape
    scape_series = [scape_series[:,i].reshape(start_scape.shape) for i in range(scape_series.shape[1])]
    #check that all the lenghts match up
    assert len(scape_series) == n, "ERROR: len(scape_series) != n"
    assert len(scape_series) == len(timesteps), "ERROR: the number of changing scapes created is not the same as the number of timesteps to be assigned to them"
    #zip the timesteps and the scapes together and return them as a list
    scape_series = list(zip(timesteps, scape_series))
    return(scape_series)


#function that takes a Landscape_Stack, a scape_num and a dictionary of {t_change:new_scape} 
#and creates an Env_Changer object that will change out Landscape_Stack.scapes[scape_num] for new_scape at
#each requisite t_change timestep
def make_custom_scape_series(scape_dict):
    pass


def make_movement_surf_series(land, scape_series, kappa=12):
    #create a deepcopy of the Landscape_Stack, to change the scape out for and then use to make movement surfaces
    copy_land = copy.deepcopy(land)
    surf_series = []
    for t, scape in scape_series:
        #change the appropriate scape in the Landscape_Stack
        copy_land.scapes[copy_land.movement_surf_scape_num].raster = scape
        #then create a Movement_Surface using the copied Landscape_Stack
        surf_series.append((t, spt.Movement_Surface(copy_land, kappa = kappa)))
    return(surf_series)


def change_scape(land, scape_num, new_scape):
    land.scapes[scape_num].raster = new_scape


def change_movement_surf(land, new_movement_surf):
    assert type(new_movement_surf) == spt.Movement_Surface, "ERROR: new_movement_surf does not appear to be of type utils.spatial.Movement_Surf"
    land.movement_surf = new_movement_surf


def get_scape_change_fn(land, scape_num, new_scape):
    def fn(land = land, scape_num = scape_num, new_scape = new_scape):
        change_scape(land = land, scape_num = scape_num, new_scape = new_scape)
    return(fn)


def get_surf_change_fn(land, scape_num, new_surf):
    def fn(land = land, scape_num = scape_num, new_surf = new_surf):
        change_movement_surf(land = land, new_movement_surf = new_surf)
    return(fn)


def get_land_change_fn(land, scape_num, item):
    type_fns = {np.ndarray              : get_scape_change_fn, 
                spt.Movement_Surface    : get_surf_change_fn
                }
    change_fn = type_fns[type(item)](land, scape_num, item)
    return(change_fn)


    ###################
    # for Pop_Changer #
    ###################
                                                           
#'bottleneck' ('instant', or 'gradual')
#'decrease'   ( rate and period of exponential decrease, or target size)
#'increase'   ( rate and period of exponential increase, or target size)
#'stochastic' ( range of fluctuations)
#'cyclical'
#'custom'

def get_dem_change_fns(pop , trajectory, rate=None, start=None, end=None, instant=False, size_target=None, size_range=None, min_size=None, max_size=None, size_format='frac'):
    if traj in ['increase', 'decrease']:
        fn = get_monotonic_dem_change_fns(pop = pop, rate = rate, start = start, end = end, instant = instant, size_target = size_target, size_format = size_format)
    elif traj == 'stochastic':
        fn = get_stochastic_dem_change_fns(instant = instant, start = start, end = end, size_range = size_range, size_format = size_format)
    elif traj == 'cyclical':
        fn = get_cyclical_dem_change_fns(start = start, end = end, size_range = size_range, min_size = min_size, max_size = max_size, size_format = size_format)
    elif traj == 'custom':
        fn = get_custom_dem_change_fns(timesteps = timesteps, sizes = sizes, size_format = size_format)
        

def make_monotonic_dem_change(pop, rate, start=None, end=None, size_target=None, size_format='frac'):
    pop.K *= rate


def get_monotonic_dem_change_fns(pop, rate, start=None, end=None, size_target=None, size_format = 'frac'):
    #get the timesteps for the demogprahic changes (subtract 1 from start to Python's 0th timestep as 1)
    #NOTE: setting start and end to the same value will create a single-timestep change (e.g. a sudden bottleneck)
    timesteps = range(start-1, end) 
    fns = []
    for t in timesteps:
        def fn(pop = pop, rate = rate, size_format = size_format):
            make_monotonic_dem_change(pop = pop, rate = rate, size_format = size_format)
        fns.append(fn)
    change_fns = list(zip(timesteps, fns))
    return(change_fns)




    ###################
    # for Gen_Changer #
    ###################


    ###################
    # for Sim_Changer #
    ###################


    ###################
    #     Others      #
    ###################
    

