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

#other imports
import numpy as np
from operator import attrgetter as ag
from collections import OrderedDict as OD
import itertools as it



#TODO:
    #1 Finish writing the get_surf_fn and get_scape_fn functions, so that they return functions ready to be
    #called to set the correct items to the correct places

    #2 review and clean the code I've written so far, and perahps restructure/simplify pieces of it

    #3 rename movement_surf to move_surf package-wide?

    #4 then think about how to further generalize these structures, and create similar ones for Gen_Changer
    #and Pop_Changer

    #5 then write the overall Changer class to unite everything

    #6 try running and debugging what I've got so far, and reflect on if it's a decent wa to do this or if it
    #needs retooling



######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

class Land_Changer:
    def __init__(self, change_scapes, params):
        #set the type-label of the changer
        self.type = 'land'

        #create self.next_change, which defaults to None but will be set by self.set_changes and upkept by self.make_change
        self.next_change = None

        #an iterator over a list (timestep, change_function) tuples, where the change function is an
        #instantiation of self.change_scape or self.change_movement_surf
        self.changes = self.set_changes(params)
       
                #create the storage dictionary, which will have keys 'scape' and 'movement_surf';
            #the 'scapes' value will be a list of (scape_num, scape) tuples (the first an integer, the second a np.array);
            #the 'movement_surfs' value will be a list of 
        self.storage = {}

    def set_next_change(self):
        self.next_change = next(self.changes)

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
                surf_series = make_movement_surf_series(scape_series)
                surfs[scape_num] = surf_series
        
        #check that timesteps are the same for the surf_series and scape_series items pertaining to the
        #movement-surface scape
        if land.movement_surf_scape_num in scape_series.keys():
            assert [i[0] for i in scape_series[land.movement_surf_scape_num]] == [i[0] for i in surf_series[land.movement_surf_scape_num]], 'ERROR: scape_series and surf_series do not contain the same timesteps for the movement surface scape'
        
        #get everything in chronological order (with scape-change functions coming before surf-change
        #functions if they're in the same timestep)
        all_changes = []
        for scape_num in scapes.keys():
            all_changes.extend([(t, scape_num, scape) for t,scape in scapes[scape_num].items()])
            if scape_num == land.movement_surf_scape_num:
                all_changes.extend([(t, scape_num, surf) for t,surf in surfs[scape_num].items()])
        all_changes = sorted(all_changes, key = lambda x: x[0])

        #make a list of change_fns for the changes in all_changes 
        change_fns = [self.get_change_fn(land, *change) for change in all_changes]
        change_fns = iter(change_fns)
        self.changes = change_fns
        self.set_next_change()

    def change_scape(self, land, scape_num, new_scape):
        land.scapes[scape_num] = new_scape

    def change_movement_surf(self, land, new_movement_surf):
        assert type(new_movement_surf) == utils.spatial.Movement_Surf, "ERROR: new_movement_surf does not appear to be of type utils.spatial.Movement_Surf"
        land.movement_surf = new_movement_surf

    def get_scape_fn(self, land, t, scape_num, scape):
        pass

    def get_surf_fn(self, land, t, scape_num, surf):
        pass

    def get_change_fn(self, land, t, scape_num, item):
        if type(item) is utils.spt.Movement_Surface:
            fn = self.get_surf_fn(t = t, scape_num = scape_num, surf = item)
        elif type(item) is np.ndarray:
            fn = self.get_scape_fn(t = t, scape_num = scape_num, scape = item)
        return(fn)

    def make_change(self, t):
        #if this is the right timestep for the next change
        if t == self.next_change[0]:
            #call the next_change function to make the change
            self.next_change[1]()
            #then load the next change after that
            self.set_next_change()
            #then recurse this function, so that multiple changes will happen this timestep if need be
            self.make_change(t)


class Pop_Changer:
    def __init__(self, change_vals, params):
        self.type = 'pop'


class Gen_Changer:
    def __init__(self, change_vals, params):
        self.type = 'gen'

#an overall Sim_Changer class, which will contain and operate the Land_, Pop_, and Gen_Changer objects
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

#function that takes a starting scape, an ending scape, a number of timesteps, and a Model object,
#and returns a linearly interpolated stack of scapes
def make_linear_scape_series(start_scape, end_scape, t_start, t_end, n):
    #get (rounded) evenly spaced timesteps at which to implement the changes
    timesteps = np.round(np.linspace(t_start, t_end, n))
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
    #zip the timesteps and the scapes together, and return them
    scape_series = zip(timesteps, scape_series)
    return(scape_series)


#function that takes a Landscape_Stack, a scape_num and a dictionary of {t_change:new_scape} 
#and creates an Env_Changer object that will change out Landscape_Stack.scapes[scape_num] for new_scape at
#each requisite t_change timestep
def make_custom_scape_series(scape_dict):
    pass



def make_movement_surf_series(land, scape_series, kappa=12):
    copy_land = copy.deepcopy(land)
    surf_series = []
    for step, scape in scape_series:
        land.scapes[land.movement_surf_scape_num] = scape
        surf_series.append((t, spt.Movement_Surface(land, kappa = kappa)))
    return(surf_series)




