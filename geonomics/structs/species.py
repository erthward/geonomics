#!/usr/bin/python
# species.py

# flake8: noqa

'''
Defines the Species class, with its associated methods and supporting functions
'''

#geonomics imports
from geonomics.utils import viz, spatial as spt
from geonomics.structs.genome import (_make_genomic_architecture,
                                      _set_empty_genomes,
                                      _make_starting_mutations)
from geonomics.structs.landscape import Layer
from geonomics.structs.individual import Individual, _make_individual
from geonomics.ops.movement import _do_movement, _do_dispersal
from geonomics.ops.mating import _find_mates, _draw_n_births, _do_mating
from geonomics.ops.selection import _calc_fitness
from geonomics.ops.mutation import _do_mutation
from geonomics.ops.demography import _do_pop_dynamics, _calc_logistic_soln
from geonomics.ops.change import _SpeciesChanger
from geonomics.sim import burnin

#other imports
import numpy as np
from numpy import random as r
import random
import matplotlib as mpl
viz._check_display()
import matplotlib.pyplot as plt
from scipy.stats.distributions import norm
from collections import Counter as C
from collections import OrderedDict as OD
import tskit
import msprime
from copy import deepcopy
from operator import itemgetter
from operator import attrgetter
import sys


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################


#a simple little class into which params-file parameters and their values will
#be dumped (as attributes). This class will be added as a hidden attribute to
#the Species. Can do this so that these parameters' values can still be 
#accessed by Geonomics functions easily (in fact, the same as if they were
#Species attributes because of how the Species' __getattr__ method is
#altered below), but they will be hidden to the regular user (so that they
#don't appear as though they could be just be changed in a custom model script
#to change model behavior (since some of them are used to instantiate objects
#at the time the model is made, so are 'baked in', such that those objects
#would need to be reinstantiated each time the params concerned were changed).
class _ParamsVals:
    def __init__(self, spp_name):
        self.spp_name = spp_name


#the Species class
class Species(OD):

    #######################
    ### SPECIAL METHODS ###
    #######################

    def __init__(self, name, idx, inds, land, spp_params,
                                                genomic_architecture=None):

        #check the inds object is correct
        assert type(inds) in (OD, dict), ("Argument inds must be of "
                                "type dict or type collections.Ordered_Dict")
        assert list(set([i.__class__.__name__ for i in list(inds.values(
            ))])) == [ 'Individual'], ("Argument inds must be a "
            "dictionary containing only instances of the individual."
                                                    "Individual class")

        #attribute to hold the Species' idx in the Community dictionary
        self.idx =  idx
        # update with the input dict of all individuals,
        #as instances of individual.Individual class
        self.update(inds)

        #set other attributes
        self.name = str(name)
        self._land_dim = land.dim
        self._land_res = land.res
        self._land_res_ratio = land._res_ratio
        self._land_ulc = land.ulc
        self._land_prj = land.prj
        #attribute to keep track of iteration number this spp is being used for
        #(optional; will be set by the iteration.py module if called)
        self._it = None
        #attribute to track total number of timesteps to be run
        #(will be set by the Model object)
        #self.T = None 
        # attribute to keep of track of number of timesteps
        #the species has evolved (starts at -1, 
        #to indicate unrun, and so that first timestep will
        #be set to 0 at beginning of timestep)
        # NOTE: This way, if model is stopped, changes 
        #are made, then it is run further,
        # this makes it easy to continue tracking from the beginning
        self.t = -1
        #will be switched to True when the species passes the burn-in tests
        self.burned = False
        #will be switched to True if the species goes extinct
        self.extinct = False
        #starting pop size
        self.start_N = len(self)
        #create a tracker for the maximum Individual index already used
        #(this will be updated each time offspring are created, so that
        #indices are never repeated, ovewrriting existing Individuals)
        self.max_ind_idx = max([*inds])
        #attribute to hold a landscape.Layer object of
        #the current species density
        self.N = None
        # attribute to hold an landscape.Layer object of
        #the local carrying capacity (i.e.
        #'target' dynamic equilibrium species density)
        self.K = None
        #the index number of the Layer to be used as this Species' K-Layer
        self.K_layer = None
        #the factor by which this Species' K-Layer should be
        #multiplied to get its K raster
        self.K_factor = None
        # list to record species size (appended each time
        self.Nt = []
        # spp.increment_age_stage() is called)
        # tracker of number of births each time spp.do_mating is called
        self.n_births = []
        # tracker of number of deaths each time
        #demography.spp_dynamics is called
        self.n_deaths = []
        #attributes for storing numpy arrays of all individuals'
        #coordinates and cells, to avoid repeat compute time each turn
        self.coords = None
        self.cells = None
        #create empty attributes to hold spatial objects that will
        #be created after the species is instantiated
        self._kd_tree = None
        #create empty attribute to hold the _DensityGridStack
        self._dens_grids = None
        #create an attribute to indicate whether this species
        #will have movement; set to False now but will be updated
        #below if a 'movement' section is encountered in the params
        self._move = False
        #create empty attributes for spatial._ConductanceSurface objects
        #that could be used for movement and/or dispersal
        #(may be created, depending on paramters)
        self._move_surf = None
        self._disp_surf = None
        #create an empty changer attribute, which will
        #be reset if the parameters define changes for this spp
        self._changer = None
        #set the sex_ratio to 0.5 default (but this will be changed if a
        #non-1/1 sex ratio is provided in the params)
        self.sex_ratio = 0.5
        #create a private _ParamsVals object as the _pv attribute
        self._pv = _ParamsVals(self.name)
        #then grab all of the mating, mortality, and movement
        #parameters as attributes of that _ParamsVals object
        for section in ['mating', 'mortality', 'movement']:
            if section in [*spp_params]:
                for att,val in spp_params[section].items():
                    #leave out the move_surf and disp_surf components,
                    #which will be handled separately
                    if not isinstance(val, dict):
                        #convert sex ratio to the probabilty of drawing a male
                        if att == 'sex_ratio':
                            val = val / (val + 1)
                        #add as an attribute of the _ParamsVals object (_pv)
                        setattr(self._pv, att, val)
                #if the section is 'movement', and it's in the params,
                #this means the Species should have movement,
                #so update the self._move attribute
                if section == 'movement':
                    if spp_params[section].move:
                        self._move = True

        #if sex is True and repro_age is an int or float, coerce to a tuple
        #(one val for each sex)
        if self.sex:
            if type(self.repro_age) in [float, int]:
                self.repro_age = (self.repro_age, self.repro_age)

        #set the GenomicArchitecture object
        self.gen_arch = genomic_architecture
        assert (self.gen_arch.__class__.__name__ == 'GenomicArchitecture'
            or self.gen_arch is None), ("The Species.gen_arch attribute "
            "must be an instance of the genome.GenomicArchitecture class "
            "or else None.")

        # create a tskit.TableCollection object, if the species uses genomes
        if self.gen_arch is not None:
            self._tc = tskit.TableCollection(sequence_length=self.gen_arch.L)
        else:
            self._tc = None

        # create an attribute to hold a TreeSequence object
        self._ts = None

        #set the selection attribute, to indicate whether or not
        #natural selection should be implemented for the species
        self.selection = (self.gen_arch is not None and
            (self.gen_arch.mu_delet > 0 or self.gen_arch.traits is not None))

        #set the self.mutate attribute (a boolean indicating whether
        #or not to enact mutation, which is True if gen_arch._mu_tot > 0
        self.mutate = (self.gen_arch is not None
                       and self.gen_arch._mu_tot is not None
                       and self.gen_arch._mu_tot > 0)

        #set the self.mut_log attribute, which dictates whether
        #or not a mutation log should be written for this spp
        self.mut_log = None
        if 'gen_arch' in [*spp_params]:
            self.mut_log = spp_params.gen_arch.mut_log

        #create a coord attrgetter function,
        #for use in getting all individs' coordinates
        self._coord_attrgetter = attrgetter('x', 'y')


    #override the __deepcopy__ method, so that the species can
    #be copy.deepcopy'd (because otherwise this doesn't work for
    #classes that inherit from #collections.OrderedDict)
    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        for k, v in self.items():
            result[deepcopy(k, memo)] = deepcopy(v, memo)
        return result


    #define the __str__ and __repr__ special methods
    #NOTE: this is not really a great representation; the Python
    #docs indicate that __repr__ should ideally
    #provide a representation could be used to recreate the
    #object, but if that is not possible then it 
    #should at least provide a string of the form '<... some
    #useful description ...>; I've attempted to do
    #the latter, inspired by a combo of what I've seen in a
    #few other packages (pandas, netCDF4, osgeo) 
    #(though perhaps I should more closely consider how to handle 
    #the params in a more nuanced, precise way once I'm
    #done writing the codebase, because right now this leaves out
    #params that have lists and dictionarires at
    #values, and also includes some parameters that are input in
    #params.py but not all, and some params that
    #are created internally but not all
    def __str__(self):
        #get a string representation of the class
        type_str = str(type(self))
        #get a string representation of the first and last individuals
        inds_str = '%i Individuals:\n\t' % len(self)
        first_ind_str = OD(self.items()).__str__().split(
                                                    '), ')[0] + ')\n\t...\n\t'
        last_ind_str = ', '.join(self.items().__str__().split(', ')[-2:])
        inds_str = inds_str + first_ind_str + last_ind_str + '\n'
        #get a string representation of the first two and last two parameters
        #params = sorted([str(k) + ': ' +str(v) for k,v in vars(
        #   self).items() if type(v) in (str, int, bool,
        #                               float)], idx = lambda x: x.lower())
        #params_str = "Parameters:\n\t" + ',\n\t'.join(params[:2]) + ','
        #params_str = params_str + '\n\t...\n\t'
        #params_str = params_str + ',\n\t'.join(params[-2:])
        #return '\n'.join([type_str, inds_str, params_str])
        return '\n'.join([type_str, inds_str])

    def __repr__(self):
        repr_str = self.__str__()
        return repr_str

    #customize the __getattr__ special method, so that attributes inside the
    #private _ParamsVals attribute (_pv) also behave as though their attributes
    #of the Species
    def __getattr__(self, attr):
        try:
            val = self._pv.__getattribute__(attr)
            return val
        except Exception:
            raise AttributeError("The Species has no attribute %s" % attr)


    #####################
    ### OTHER METHODS ###
    #####################

        #################
        #private methods#
        #################

    #method to calculate and set self.K
    def _set_K(self, land):
        self.K = land[self.K_layer].rast * self.K_factor

    #method to set self.N
    def _set_N(self, N):  #NOTE: Requires a landscape.Layer instance
        self.N = N

    #method to append current spp size to the spp.Nt list
    def _set_Nt(self):
        self.Nt.append(len(self))

    #method to increment the self.t attribute (i.e. the timestep counter)
    def _set_t(self):
        self.t += 1

    #method to reset the self.t attribute (i.e. the timestep counter)
    def _reset_t(self):
        self.t = -1

    #method to increment all species' age by one
    #(also adds current spp size to tracking array)
    def _set_age_stage(self):
        # increment age of all individuals
        [ind._set_age_stage() for ind in self.values()];

    #method to move all individuals simultaneously, and sample
    #their new locations' environment
    def _do_movement(self, land):
        _do_movement(self)
        self._set_e(land)
        self._set_coords_and_cells()

    #function for finding all the mating pairs in a species
    def _find_mating_pairs(self):
        mating_pairs = _find_mates(self, sex=self.sex,
                                   dist_weighted_birth=self.dist_weighted_birth)
        return mating_pairs

    #function for executing mating for a species
    def _do_mating(self, land, mating_pairs, burn=False):

        #draw the number of births for each pair, and append
        #total births to self.n_births list
        if self.n_births_fixed:
            n_births = np.array(
                            [self.n_births_distr_lambda] * len(mating_pairs))
        else:
            n_births = _draw_n_births(len(
                                    mating_pairs), self.n_births_distr_lambda)
        total_births = sum(n_births)
        self.n_births.append(total_births)

        #create the offspring_ids
        next_offspring_key = self.max_ind_idx + 1
        offspring_keys = set(range(next_offspring_key,
                                        next_offspring_key + total_births))
        #update self.max_ind_idx
        if len(offspring_keys) > 0:
            self.max_ind_idx = max(offspring_keys)

        #copy the keys, for use below in do_mutation, _set_z, _calc_fitness
        keys_list = [*offspring_keys]

        if not burn and self.gen_arch is not None:
            recomb_events = self.gen_arch._recomb_events._get_events(
                                                            total_births*2)
            # NOTE: this gives us a hierarchically organized list,
            # containing both the offspring's new genomes and their
            # segment information (to be added to the
            # tskit.TableCollection.edges table), organized as a sublist for
            # each pair, then a subtuple within that, containing a new genome
            # and then a 2-tuple of segments info for each of the new
            # offspring's homologues
            seg_sets = _do_mating(self, mating_pairs, n_births, recomb_events)

        for n_pair, pair in enumerate(mating_pairs):

            parent_midpoint_x = (self[pair[0]].x + self[pair[1]].x)/2
            parent_midpoint_y = (self[pair[0]].y + self[pair[1]].y)/2

            n_offspring = n_births[n_pair]

            for n in range(n_offspring):

                #get the next offspring_key
                offspring_key = offspring_keys.pop()

                offspring_x, offspring_y = _do_dispersal(
                    self, parent_midpoint_x, parent_midpoint_y,
                    self.dispersal_distr_mu, self.dispersal_distr_sigma)

                #set the age to 0
                age = 0

                #set the sex correctly
                if self.sex:
                    sex = r.binomial(1, self.sex_ratio)
                else:
                    sex = None

                #set the new_genome correctly
                if self.gen_arch is not None:
                    if burn:
                        pass
                        #new_genome = np.array([0])
                    else:
                        #new_genome = genomes_and_segs[n_pair][n][0]
                        segs = seg_sets[n_pair][n]
                else:
                    #new_genome = None
                    pass

                #create the new individual
                self[offspring_key] = Individual(idx = offspring_key,
                                                 age = age, x = offspring_x,
                                                 y = offspring_y, sex = sex)

                # during the main phase, for species with genomes,
                # update the tskit tables
                if (self.gen_arch is not None
                    and not burn):
                    loc = [offspring_x, offspring_y]
                    #if self.gen_arch.traits is not None:
                    #    loc = (loc + [*self[offspring_key].z] +
                    #        [self[offspring_key].fit])
                    #print(loc)
                    #print(self[offspring_key].fit)
                    offspring_ind_id = self._tc.individuals.add_row(
                        location=loc,
                        # NOTE: using the metadata column to store to gnx
                        # individual idx, for later matching to update
                        # Individual._individuals_tab_id after tskit's simplify
                        # algorithm filters individuals
                        metadata=offspring_key.to_bytes(length=4,
                                                        byteorder='little'))
                    self[offspring_key]._individuals_tab_id = offspring_ind_id

                    # add rows to the nodes table, setting
                    # the 'flags' column vals to 0
                    # (to indicate they're not considered sample nodes),
                    # and setting the 'individual' column vals to ids
                    # returned from tc.individuals.add_row(), then adding
                    # the returned tskit Node ids to Individual_nodes_tab_ids
                    # attribute (which is a list)
                    # NOTE: make time negative so that parent time is always
                    # greater than child time (as it would be expressed in the
                    # coalescent, except that we can't use positive numbers
                    # here because we want to allow for the possibility
                    # that a model could be walked for any arbitrary
                    # number of time steps)
                    self[offspring_key]._set_nodes_tab_ids(
                        *[self._tc.nodes.add_row(flags=1, time=-self.t,
                                                 population=0,
                        individual=offspring_ind_id) for _ in range(
                                                            self.gen_arch.x)])

                    # add edges to the tskit edges table
                    # NOTE: `[*seg_set]` is necessary to star-unpack
                    # the zip object
                    edge_id = [self._tc.edges.add_row(
                        parent=seg[0], left=seg[1], right=seg[2],
                        child=self[offspring_key]._nodes_tab_ids[homol]
                        ) for homol, seg_set in enumerate(
                        segs) for seg in [*seg_set]]

        # do mutation if necessary
        if self.mutate and not burn:
             _do_mutation(keys_list, self, log = self.mut_log)

        # sample all individuals' environment values, to initiate for offspring
        self._set_e(land)
        self._set_coords_and_cells()

        #set new individuals' phenotypes (won't be set
        #during burn-in, because no genomes assigned;
        #won't be set if the species has no gen_arch)
        if (self.gen_arch is not None and self.gen_arch.traits is not None
            and not burn):

            # update the TreeSequence object
            self._tc.sort()
            self._ts = self._tc.tree_sequence()
            #vars = self._tc.tree_sequence().variants(samples=

            self._set_z(individs=keys_list)
            self._calc_fitness(set_fit=True, individs=keys_list)


    #method to do species dynamics
    def _do_pop_dynamics(self, land):
        #implement selection, iff self.selection is True and the spp has
        #already been burned in
        with_selection = self.selection and self.burned
        burn = not self.burned
        #then carry out the pop-dynamics, with selection as set above, and save
        #result, which will be True iff spp has gone extinct
        extinct = _do_pop_dynamics(self, land,
            with_selection = with_selection, burn = burn)
        if extinct:
            #set self.extinct equal to True, so that the iteration will end
            self.extinct = extinct

    #method to make species changes
    def _make_change(self, verbose=False):
        self._changer._make_change(t = self.t, additional_args = {'spp': self},
            verbose = verbose)

    #method to check if the species has gone extinct
    def _check_extinct(self):
        return len(self) == 0

    #method to calculate species density
    def _calc_density(self, normalize = False, as_layer = False, set_N=False):

        '''
        Calculate an interpolated raster of local species density, using
        the spatial._DensityGridStack object stored in the
        Species._dens_grids attribute.

        If normalize is True, the density raster will vary between 0 and 1
        (after being normalized by its own maximum value). If false, it will
        vary between 0 and its maximum estimated density value.
        '''
        #check validity of normalize argument
        assert type(normalize) is bool, ("The 'normalize' argument takes "
            "a boolean value.\n")
        #get species coordinates
        x = self._get_x()
        y = self._get_y()
        #calculate the density array
        dens = self._dens_grids._calc_density(x, y)
        #set min value to 0
        dens = np.clip(dens, a_min = 0, a_max = None)
        #normalize, if necessary 
        if normalize:
            # Use max_val to normalize the density raster to either
            #0 to its current max val or 0 to 1, to make sure 
            #the interpolation didn't generate any values slightly
            #outside this range
            norm_factor = dens.max() - dens.min()
            dens = (dens - dens.min()) / norm_factor
        #return as layer, if necessary
        if as_layer == True:
            dens = Layer(dens, 'density', 'density', self._land_dim)
        #set self.N attribute, if necessary
        if set_N:
            self._set_N(dens)
        #or else return the density array
        else:
            return dens

    #method to set the individuals' environment values 
    def _set_e(self, land, individs = None):
        if individs is None:
            inds_to_set = self.values()
        else:
            ig = itemgetter(*individs)
            inds_to_set = ig(self)
            if isinstance(inds_to_set, individual.Individual):
                inds_to_set = (inds_to_set,)
        hab = [ind._set_e([lyr.rast[int(ind.y), int(
            ind.x)] for lyr in land.values()]) for ind in inds_to_set]

    #method to set the individuals' phenotype attributes 
    def _set_z(self, individs=None):
        if individs == None:
            individs = [*self]
        z_dicts = [self._get_genotypes_or_phenotypes(
                                    with_dominance=self.gen_arch._use_dom,
                                    individs=individs, trait=trt.idx,
                                    ) for trt in self.gen_arch.traits.values()]
        [setattr(self[ind], 'z', np.hstack(
                    [d[ind] for d in z_dicts])) for ind in individs];

    #method to set the individuals' fitness attributes
    def _set_fit(self, fit, individs=None):
        if individs is None:
            individs = [*self]
        [self[ind]._set_fit(f) for ind, f in zip(individs, fit)];


    #method to set species' coords and cells arrays
    def _set_coords_and_cells(self):
        self.coords = self._get_coords()
        self.cells = np.int32(np.floor(self.coords))


    #method to set the species' kd_tree attribute (a spatial._KDTree)
    def _set_kd_tree(self, leafsize = 100):
        self._kd_tree = spt._KDTree(coords = self.coords, leafsize = leafsize)


    #method to set the species' spatial._DensityGridStack attribute
    def _set_dens_grids(self, land, widow_width = None):
        self._dens_grids = spt._DensityGridStack(land = land,
                                window_width = self.density_grid_window_width)


    # method to fill the tskit.TableCollection's tables (to be called
    # after the model has burned in
    def _set_genomes_and_tables(self, burn_T, T):
        # set the species' neutral and non-netural loci,
        # and set genomes to all zeros
        _set_empty_genomes(self, burn_T, T)

        # simulate a coalescent ancestry with number of samples equal to our
        # species' number of haploid genomes (i.e. 2*N_0)
        ts = msprime.simulate(len(self) * 2, Ne=1000, length=self.gen_arch.L)

        # then grab the simulation's tree, and the tree's TableCollection
        # NOTE: the TreeSequence object only has one tree,
        # because no recombination was used in the sim
        tree = ts.first()
        tables = ts.dump_tables()

        # set the sequence length
        #tables.sequence_length = self.gen_arch.L
        # clear the mutations table
        tables.mutations.clear()

        # loop over all sites, so that each site's row goes into
        # the sites table in order (such that 1.) there's no need
        # to track how gnx sites map onto sites-table row ids,
        # and 2.) there will be no need to deduplicate sites later on)
        for site in range(self.gen_arch.L):
            #determine whether this is a neutral or non-neutral site
            if site in self.gen_arch.neut_loci:
                metadata='n'.encode('ascii')
            elif site in self.gen_arch.nonneut_loci:
                metadata='t'.encode('ascii')
            # add the variant's site to the sites table
            tables.sites.add_row(position=site, ancestral_state='0',
                                 metadata=metadata)

        # grab the nodes flags, which are 1 for current nodes,
        # 0 for past nodes, into two separate objects
        current_nodes = [*np.where(tables.nodes.flags == 1)[0]]
        # reverse, so that I can pop nodes off the 'front'
        current_nodes = current_nodes[::-1]
        past_nodes = [*np.where(tables.nodes.flags != 1)[0]]
        # create an empty list, to fill up the individual ids for each node
        # NOTE: setting to a vector of -1 initially, to easily check that
        # all values have been assigned at the end by asking if all >= 0 
        nodes_tab_individual_col = np.int32(np.ones(len(
                                                    tables.nodes.flags))*-1)

        # NOTE: there are no requirements or restrictions for the individuals
        # and nodes tables (e.g. order, etc.), and thus the tables' order
        # is not affected by the TableCollection.simplify algorithm.
        # So, I could add either current or past individuals and nodes first;
        # choosing to add past first, so that all 'real' individuals
        # from the start of the geonomics simulation forward will
        # wind up having their individuals and nodes rows in a single
        # block at the tables' bottoms

        # add an individual to the individuals table for
        # each coalescent-simulated node before the current time
        # NOTE: adding no metadata, and no location, to indicate that this is
        # a 'fake' individual, invented just to match up to the nodes
        # simulated for the starting population
        for node in past_nodes:
            ind_id = tables.individuals.add_row(flags=0)
            # store its individual id in the nodes table's individuals column
            nodes_tab_individual_col[node] = ind_id

        # create and add to the individuals table a new row for
        # each real individual
        for ind in self.values():
            # get the 'location' column info, which will include the x and y
            # positions of an individual, as well as the individual's
            # phenotypes and fitness, if traits are being used
            loc = [ind.x, ind.y]
            if self.gen_arch.traits is not None:
                loc = loc + ind.z + [ind.fit]
            # add a new row to the individuals table, setting the location
            # column's value to loc
            # NOTE: using the metadata column to store to the gnx
            # individual idx, for later matching to update
            # Individual._individuals_tab_id after tskit's simplify
            # algorithm filters individuals
            ind_id = tables.individuals.add_row(flags=1, location=loc,
                metadata=ind.idx.to_bytes(length=4, byteorder='little'))
            ind._individuals_tab_id = ind_id

            # assign the individual 2 randomly chosen nodes from
            # the current time step, and associate the individual's
            # _individuals_tab_id with those 2 nodes in some data
            # structure to collect this
            ind_node_ids =[current_nodes.pop() for _ in range(self.gen_arch.x)]
            ind._set_nodes_tab_ids(*ind_node_ids)
            # add this individual's ind_id to the
            # nodes_tab_individual_col list, once for each node
            nodes_tab_individual_col[ind_node_ids] = ind_id
        # make sure that all nodes were assigned to individuals
        assert np.all(nodes_tab_individual_col >= 0), ('Some nodes not '
                                                       'given individs')

        nodes_cols = tables.nodes.asdict()
        # use the node_tab_individual_col data structure to reset
        # the individuals column in the nodes table
        nodes_cols['individual'][:] = nodes_tab_individual_col
        # increment all the birth times by 1, so that all individuals are
        # marked as having been born before the start of the model's main phase
        nodes_cols['time'] += 1
        tables.nodes.set_columns(**nodes_cols)

        # add sufficient mutations, at only current nodes, to produce
        # the starting 1-allele frequencies parameterized for this species
        _make_starting_mutations(self, tables)

        # TODO: ADD PROVENANCES ROW!

        # assign as the species' TableCollection
        self._tc = tables
        self._tc.sort()
        self._ts = self._tc.tree_sequence()
        
        # set individuals' phenotypes and fitness, if necessary
        if self.gen_arch.traits is not None:
            self._set_z()
            self._calc_fitness(set_fit=True)

        return


    # method to sort and simplify the TableCollection,
    # and update Individuals' node IDs
    # NOTE: rather than deduplicate sites (and sometimes
    # computer parent sites), as recommended by tskit's docs,
    # we only have to sort and simplify here, because we have
    # chosen to just add all simulated sites to the sites table (in
    # Species._set_table_collection) at the simulation's outset,
    # and also to only allow an infinite-sites model
    # (FOR NOW, anyhow... could relax these constraints later, if desired,
    # in which case would need to revamp this approach)
    def _do_sort_simp_table_collection(self, check_nodes=False,
                                       check_individuals=False,
                                       check_genotypes=False):
        # sort the TableCollection
        self._tc.sort()
        # get an array of all the current individuals' nodes,
        # the nodes for which the tables will be simplified
        curr_nodes = np.hstack([[*ind._nodes_tab_ids.values(
                                                )] for ind in self.values()])
        # run code necessary for checing that individuals' table ids are
        # correctly assigned, if requested
        if check_individuals:
        # NOTE: this check will only work for models that have at least one
        # trait, because only those models will have at least 3 values in the 
        # location column's value (and only the values after the first 2 will
        # be the same across an individual's lifetime, since they move position
        # frequently but do not cahnge phenotype at all)

            meta_b4 = self._tc.individuals.metadata
            loc = self._tc.individuals.location
            locoff = self._tc.individuals.location_offset
            loc_b4_by_table_ids = {
                table_id: (loc[locoff[table_id]: locoff[
                    table_id+1]]) for table_id in range(len(
                                                self._tc.individuals.flags))}
            z0_b4 = {ind.idx: loc_b4_by_table_ids[
                        ind._individuals_tab_id][2] for ind in self.values()}

        # now simplify the tables and get the new ids output 
        # from the tc.simplify method
        # NOTE: we are not using the ouput, but it is an array where
        # each value is the new node ID of the node that was in that index's
        # row position in the old nodes table (and -1 if that node was
        # dropped during simplication); see tskit docs for details
        output = self._tc.simplify(curr_nodes, filter_individuals=True,
                                   filter_sites=False)

        # make an Nx3 np array containing 1.) gnx ids, 2.) the new
        # homologue 0 node ids, and 3.) the new homologue 1 node ids,
        # in each its 3 cols; then loop over its rows to update
        # each individual's _nodes_tab_ids attribute
        inds_gnx_ids = np.array([*self]).reshape((len(self), 1))
        new_ids = np.hstack((inds_gnx_ids, np.array([*range(len(
            inds_gnx_ids) * 2)]).reshape((len(inds_gnx_ids), 2))))
        for id, n0, n1 in new_ids:
            self[id]._nodes_tab_ids.update({0:n0, 1:n1})

        # update Individuals' table ids
        # (i.e. their Individual._individuals_tab_id attributes)
        ind_meta = self._tc.individuals.metadata
        #if check_individuals:
            #print('b4', len(meta_b4), 'af', len(ind_meta))
        ind_off = self._tc.individuals.metadata_offset
        new_individuals_tab_id = {int.from_bytes(
            ind_meta[ind_off[i]: ind_off[i+1]].tobytes(),
            'little'): i for i in range(len(ind_off) - 1)}
        [setattr(ind, '_individuals_tab_id',
                 new_individuals_tab_id[id]) for id, ind in self.items()]

        # check that individuals' nodes-table ids were correclty updated,
        # if the check is requested
        if check_nodes:
            import pandas as pd
            #create another identically structured Nx3 np array, to hold the
            #individuals' gnx ids and their node ids according to the tskit
            #nodes table (for cross-checking)
            new_ids_from_tables = np.ones((len(inds_gnx_ids), 3)) * np.nan
            new_ids_from_tables[:,0] = inds_gnx_ids[:,0]
            nodedf = pd.DataFrame({k:v for k, v in self._tc.nodes.asdict(
                ).items() if k in ['time', 'individual']})
            ind_table_ids = [self[idx[
                            0]]._individuals_tab_id for idx in inds_gnx_ids]
            # loop over individuals' gnx ids and corresponding table row ids,
            # filling up the comparison table
            for ind, idx in zip(inds_gnx_ids[:,0], ind_table_ids):
                subdf = nodedf[nodedf['individual'] == idx]
                #print('...\ngnx id:', ind, '\ttable id:', idx)
                #print('gnx node ids:  ', [*self[ind]._nodes_tab_ids.values()])
                #print('table node ids:', [*subdf.index])
                node_id_vals = [*subdf.index]
                row_idx = np.where(new_ids_from_tables[:,0] == ind)[0][0]
                new_ids_from_tables[row_idx, 1:] = node_id_vals
            # now make sure tabls are identical
            assert np.all(new_ids_from_tables == new_ids), ("All the node "
                                                            "ids don't match!")
            print('PASS: all node IDs appear correctly reassigned\n')

        # check that individuals' individuals-table ids were correclty updated,
        # if the check is requested
        if check_individuals:
            # check that individuals' table ids were correctly reassigned, by
            # way of assuring that individuals, as indexed by both their old
            # and their new individuals-table ids, have identical phenotypes
            loc = self._tc.individuals.location
            locoff = self._tc.individuals.location_offset
            loc_af_by_table_ids = {
                table_id: (loc[locoff[table_id]: locoff[
                    table_id+1]]) for table_id in range(len(
                                                self._tc.individuals.flags))}
            z0_af = {ind.idx: loc_af_by_table_ids[
                        ind._individuals_tab_id][2] for ind in self.values()}
            #print(z0_af)
            z0_check = [z0 == z0_b4[idx] for idx, z0 in z0_af.items()]
            for idx in [*self]:
                print(idx, ' b4 ', z0_b4[idx])
                print(idx, ' af ', z0_af[idx])
                print('..................')
            assert np.all(z0_check), ('phenotypes are not the same! '
                                      ' (MAKE SURE YOU USED A MODEL WITH'
                                      ' AT LEAST ONE TRAIT!)')
            print("PASS: all individuals' IDs appear correctly reassigned\n")
        if check_genotypes:
            whether_genotypes_match = []
            samples = [[*ind._nodes_tab_ids.values()] for ind in self.values()]
            samples = np.int32(np.hstack(samples))
            vars = self._tc.tree_sequence().variants(samples=samples)
            for site in range(self.gen_arch.L):
                gnx_genotypes = np.hstack([ind.g[site,
                                                 :] for ind in self.values()])
                tskit_genotypes = next(vars).genotypes
                res = np.all(gnx_genotypes == tskit_genotypes)
                whether_genotypes_match.append(res)
            final_res = np.all(whether_genotypes_match)
            assert final_res, "genotypes did not all match!"
            print("PASS: all individuals' genotypes appear correct\n")


    # method to get individs' environment values
    def _get_e(self, lyr_num=None, individs=None):
        if individs is None:
            if lyr_num is None:
                e = np.array([ind.e for ind in self.values()])
            else:
                e = np.array([ind.e[lyr_num] for ind in self.values()])
        else:
            ig = itemgetter(*individs)
            if lyr_num is None:
                e = {i:ind.e for i, ind in self.items()}
                e = np.array(ig(e))
            else:
                e = {i:ind.e[lyr_num] for i, ind in self.items()}
                e = np.array(ig(e))
        return e


    def _get_genotypes_or_phenotypes(self, loci=None, biallelic=False,
                                     individs=None, with_dominance=True,
                                     trait=None, as_array=False):
        # ensure that biallelic is True, if phenotypes are requested
        if trait is not None:
            biallelic = False
        # get the list of loci
        if loci is None:
            if trait is None:
                loci = [*range(self.gen_arch.L)]
            else:
                loci = self.gen_arch.traits[trait].loci
        # get the list of individuals
        if individs is None:
            individs = [*self]
        # and get the list of their nodes
        samples = np.int32(np.hstack([[*self[ind]._nodes_tab_ids.values(
                                                    )] for ind in individs]))
        assert len(samples) == self.gen_arch.x * len(individs), ('Number of '
                                                                 'nodes does '
                                                                 'not match '
                                                                 'number of '
                                                                 'individs!')

        # sort the TableCollection and get the TreeSequence
        #self._tc.sort()
        #ts = self._tc.tree_sequence()

        # get the genotypes for requested individuals at all requested loci
        genotypes = [var.genotypes for n, var in enumerate(self._ts.variants(
                                                samples=samples)) if n in loci]
        genotypes = np.vstack(genotypes)
        genotypes = {ind: genotypes[:, i: i + 2] for ind, i in zip(
                individs, [*range(0, len(self) * self.gen_arch.x, 2)])}

        # calculate non-biallelic or dominance-based genotypes, if needed
        if not biallelic:
            genotypes = {ind: np.mean(gt,
                                      axis=1) for ind, gt in genotypes.items()}
            if with_dominance == True:
                d = self.gen_arch.dom[loci]
                genotypes = {ind: np.clip(gt * (1 + d), a_min=None,
                                   a_max=1) for ind, gt in genotypes.items()}
        # multiply by an array of factors, if requested (i.e. in order to
        # calculate effect sizes, for phenotype calculation)
        if trait is not None:
            if len(self.gen_arch.traits[trait].loci) > 1:
                phenotypes = {ind: 0.5 + np.sum(
                         gt * self.gen_arch.traits[trait].alpha) for ind,
                         gt in genotypes.items()}
            else:
                phenotypes = {ind: gt[0] for ind, gt in genotypes.items()}
            # return as an N x L array (or 2N x L, if biallelic),
            # rather than a dict, if requested
            if as_array:
                if biallelic:
                    phenotypes = np.hstack([pt for pt in phenotypes.values(
                                                                        )]).T
                else:
                    phenotypes = np.vstack([pt for pt in phenotypes.values()])
            return phenotypes
        # return as an N x L array (or 2N x L, if biallelic),
        # rather than a dict, if requested
        if as_array:
            if biallelic:
                genotypes = np.hstack([gt for gt in genotypes.values()]).T
            else:
                genotypes = np.vstack([gt for gt in genotypes.values()])
        return genotypes

    #convenience method for getting a scalar attribute for some or all individs
    def _get_scalar_attr(self, attr_name, individs=None):
        if individs is None:
            vals = np.array([getattr(ind, attr_name) for ind in self.values()])
        else:
            ig = itemgetter(*individs)
            vals = {i: getattr(ind, attr_name) for i, ind in self.items()}
            vals = np.array(ig(vals))
        return vals

    #convenience method for getting age of whole species
    def _get_age(self, individs=None):
        ages = self._get_scalar_attr('age', individs=individs)
        return ages

    # convenience method for getting whole species' phenotype
    def _get_z(self, individs=None):
        zs = self._get_scalar_attr('z', individs=individs)
        return zs

    #convenience method for getting whole species' fitnesses
    def _get_fit(self, individs = None):
        fits = self._get_scalar_attr('fit', individs=individs)
        return fits

    def _calc_fitness(self, trait_num=None, set_fit=True, individs=None):
        fit = _calc_fitness(self, trait_num=trait_num, individs=individs)
        #set individuals' fitness attributes, if indicated
        if set_fit:
            self._set_fit(fit, individs=individs)
        return fit

    def _get_dom(self, locus):
        return {locus: self.gen_arch.h[locus]}

    def _get_coords(self, individs=None, as_float=True):
        coords = list(map(self._coord_attrgetter, self.values()))
        if individs is not None:
            ig = itemgetter(*individs)
            coords = ig(dict(zip([*self], coords)))
        if as_float:
            coords = np.float64(coords)
        else:
            coords = np.int32(np.floor(coords))
        # make sure it's at least 2d (in case a single individual is requested)
        coords = np.atleast_2d(coords)
        return coords

    def _get_plot_coords(self, individs=None, cell_coords=False):
        coords = self._get_coords(individs=individs)
        if not cell_coords:
            coords[:, 0] = coords[:, 0] * self._land_res[0] + self._land_ulc[0]
            coords[:, 1] = coords[:, 1] * self._land_res[1] + self._land_ulc[1]
        return coords

    def _get_cells(self, individs=None):
        cells = self._get_coords(individs=individs, as_float=False)
        return cells

    def _get_x(self, individs=None):
        coords = self._get_coords(individs=individs)
        return coords[:, 0]

    def _get_y(self, individs=None):
        coords = self._get_coords(individs=individs)
        return coords[:, 1]

    # method to return an n-length list of random individs;
    # return individuals, or indices, as indicated
    def _get_random_individuals(self, n, return_format='index'):
        assert return_format in [
            'index', 'individual'], ("Argument return_format can take only "
                                     "'index' or 'individual' as values "
                                     "(defaults to 'index').")
        choices = choices = r.choice(list(range(len(self))), n)
        inds = np.array(list(self.keys()))[choices]
        if return_format=='individual':
            inds = [self[ind] for ind in inds]
        return(inds)

    # method to reduce a spp to some chosen number, by randomly removing
    # N_curr_t - n individuals
    def _reduce(self, n):
        inds = [*self]
        keep = np.random.choice(inds, n, replace=False)
        for ind in inds:
            if ind not in keep:
                self.pop(ind)

    #use the kd_tree attribute to find nearest neighbors either
    #within the species, if within == True, or between the species
    #and the points provided, if within == False and points is not None
    def _find_neighbors(self, dist, within=True, coords=None, k = 2):
        #NOTE: In lieu of a more sophisticated way of
        #determining whether the kd_tree needs to be updated 
        #(i.e. if the species have undergone movement, mating,
        #or mortality since the last time it was 
        #constructed), and in an effort to minimize the number
        #of times it is constructed each time (since
        #it's not an inconsequential runtime, so telling the
        #model to rebuild it after each move, birth, or
        #death step could be unncessarily costly), for now I
        #am just telling the tree to be rebuilt each time 
        #the spp.find_neighbors() method is called!
        self._set_kd_tree()
        #if neighbors are to be found within the species,
        #set coords to self.coords (otherwise, the coords to
        #find nearest neighbors with should have been provided)
        if within:
            coords = self.coords
        #query the tree
        dists, pairs = self._kd_tree._query(coords = coords,
                                                        dist = dist, k = k)
        return(dists, pairs)

    # method for plotting the species (or a subset of its individuals, by ID)
    # on top of a layer (or landscape)
    def _plot(self, lyr_num=None, land=None, hide_land=False, individs=None,
              text=False, color='black', edge_color='face', text_color='black',
              cbar=True, size=25, text_size=9, land_cmap=None, pt_cmap=None,
              alpha=False, zoom_width=None, x=None, y=None, vmin=None,
              vmax=None, ticks=None, mask_rast=None, animate=False,
              cell_coords=False):
        # convert individs to a list (in case comes in as a numpy array)
        if individs is not None and not isinstance(individs, list):
            individs = list(individs)
        # get coords
        coords = self._get_plot_coords(individs=individs,
                                       cell_coords=cell_coords)
        # get text
        if text:
            if individs is None:
                text = [*self]
            else:
                text = individs
        else:
            text = None
        # set the plt_lims
        plt_lims = viz._get_plt_lims(land, x, y, zoom_width)
        # plot the layer(s)
        if hide_land:
            pass
        else:
            # get the layers' vmin and vmax values, if any of the layers to
            # be plotted has a change event
            if ((lyr_num is None and land._changer is not None) or
                (land._changer is not None
                 and lyr_num in [*land._changer.change_info])):
                if lyr_num is None:
                    land_vmin = [lyr._scale_min for lyr in land.values()]
                    land_vmax = [lyr._scale_max for lyr in land.values()]
                else:
                    land_vmin = [land[lyr_num]._scale_min]
                    land_vmax = [land[lyr_num]._scale_max]
            else:
                land_vmin = land_vmax = None
            viz._plot_rasters(land, lyr_num=lyr_num,
                              cbar=cbar, cmap=land_cmap, plt_lims=plt_lims,
                              ticks=ticks, mask_rast=mask_rast,
                              vmin=land_vmin, vmax=land_vmax)
        # and plot the individuals
        points = viz._plot_points(coords, lyr_num=lyr_num, color=color,
                                  edge_color=edge_color,
                                  text_color=text_color, size=size,
                                  text_size=text_size, alpha=alpha,
                                  text=text, plt_lims=plt_lims,
                                  pt_cmap=pt_cmap, vmin=vmin, vmax=vmax,
                                  animate=animate)
        return points


    #method for plotting the species on top of its estimated
    #species-density raster
    def _plot_density(self, land, normalize=False, individs=None,
            text=False, color='black', edge_color='face',
            text_color='black', size=25, text_size = 9,
            alpha=0.5, zoom_width=None, x=None, y=None, ticks=None,
            mask_rast=None):
        assert type(normalize) is bool, ("The 'normalize' argument takes "
            "a boolean value.\n")
        #update the species' coordinates and cells, in case it hasn't
        #been update since some internal or manual changes in population-size
        #have occurred
        self._set_coords_and_cells()
        dens = self._calc_density(normalize = normalize)
        plt_lims = viz._get_plt_lims(land, x, y, zoom_width)
        if normalize:
            viz._plot_rasters(dens, plt_lims = plt_lims, lyr_name = 'density',
                              ticks=ticks, mask_rast=mask_rast)
        else:
            viz._plot_rasters(dens, plt_lims = plt_lims, vmax = dens.max(),
                lyr_name = 'density', ticks=ticks, mask_rast=mask_rast)
        self._plot(hide_land=True, individs = individs, text = text,
            color=color, edge_color = edge_color, text_color = text_color,
            size=size, text_size = text_size, alpha=alpha,
                                    zoom_width = zoom_width, x = x, y = y)


    # method for plotting individuals colored by their genotype at a locus
    def _plot_genotype(self, locus, lyr_num=None, individs=None,
                       text=False, size=25, text_size = 9, edge_color='black',
                       text_color='black', cbar=True, alpha=1,
                       with_dominance=False, zoom_width=None, x=None, y=None,
                       ticks=None, mask_rast=None):

        if with_dominance == True:
            genotypes = self._get_genotypes_or_phenotypes([locus],
                                                          with_dominance=True)
        else:
            genotypes = self._get_genotypes_or_phenotypes([locus])

        if individs is not None:
            genotypes = {i:v for i,v in genotypes.items() if i in individs}

        # just assign black, gray, and white (since there's no reason
        # necessarily that these should be mapped to a certain layer, the way
        # phenotype should
        colors = ['#000000', '#808080', '#FFFFFF']

        for n, genotype in enumerate([0.0, 0.5, 1.0]):
            genotype_individs = [i for i, g in genotypes.items(
                                        ) if np.atleast_1d(g)[0] == genotype]
            # plot if there are any individuals of this genotype
            if len(genotype_individs) >= 1:
                self._plot(lyr_num = lyr_num, individs = genotype_individs,
                    text = text, color = colors[n], edge_color = edge_color,
                    text_color = text_color, cbar = cbar,
                    size = size, text_size = text_size, alpha = alpha,
                    zoom_width = zoom_width, x = x, y = y, vmin = 0, vmax = 1,
                    ticks=ticks, mask_rast=mask_rast)


    # method for plotting individuals colored by their phenotypes
    #for a given trait
    def _plot_phenotype(self, trait, lyr_num=None, land=None,
                        individs=None, text=False, size=25, text_size=9,
                        edge_color='black', text_color='black', cbar=True,
                        alpha=1, zoom_width=None, x=None, y=None,
                        ticks=None, mask_rast=None, animate=False):

        # get the trait's lyr_num, if no lyr_num provided
        lyr_num = self.gen_arch.traits[trait].lyr_num

        #z = OD(zip([*self], self._get_z()[:, trait]))
        z = self._get_genotypes_or_phenotypes(trait=trait, individs=individs,
                                              as_array=True)
        #if individs is not None:
        #    z = {i: v for i, v in z.items() if i in individs}

        # get the correct cmap for this trait's layer
        pt_cmap = viz._choose_cmap(self.gen_arch.traits[trait].lyr_num)

        points = self._plot(lyr_num = lyr_num, land = land,
                            individs = individs, text = text,
                            color = z, pt_cmap = pt_cmap,
                            edge_color = edge_color, text_color = text_color,
                            cbar = cbar, size = size, text_size = text_size,
                            alpha = alpha, zoom_width = zoom_width, x = x,
                            y = y, vmin = 0, vmax = 1, ticks=ticks,
                            mask_rast=mask_rast, animate=animate)

        return points


    # method for plotting individuals colored by their overall fitnesses,
    #or by their fitnesses for a single trait (if trait is not None)
    def _plot_fitness(self, trt_num=None, lyr_num=None, land=None,
                      individs=None, text=False, phenotype_text=False,
                      phenotype_text_color='black', fitness_text=False,
                      fitness_text_color='#333333', size=100, text_size = 9,
                      edge_color='black', text_color='black',
                      fit_cmap = 'RdYlGn', cbar=True, fitness_cbar=True,
                      alpha=1, zoom_width=None, x=None, y=None, ticks=None,
                      mask_rast=None):

        #return messages if species does not have genomes or traits
        if self.gen_arch is None:
            print(("Species._plot_fitness is not valid for species "
                   "without genomes.\n"))
            return
        elif self.gen_arch.traits is None:
            print(("Species._plot_fitness is not valid for species "
                   "without traits.\n"))
            return

        # get the trait's lyr_num, if lyr_num wasn't provided but trt_num was
        if trt_num is not None and lyr_num is None:
            lyr_num = self.gen_arch.traits[trt_num].lyr_num

        # get all individs' fitness values,
        # and get appropriate colormap
        if trt_num is None:
            w = self._calc_fitness()
            pt_cmap = 'Greys_r'
        else:
            w = self._calc_fitness(trait_num = trt_num)

        #filter out unwanted individs, if necessary
        w = OD(zip([*self], w))
        if individs is not None:
            w = {i:v for i,v in w.items() if i in individs}

        # calc minimum possible fitness (for phenotypes within 0 <= z <= 1,
        #which in reality isn't a constraint, but values lower than
        #this will also be constrained to the minimum-value color for plotting)
        #NOTE: the np.atleast_2d(...).min() construct makes this
        #work both for fixed and spatially varying phi
        if trt_num is None:
            min_fit = np.product([1 - np.atleast_2d(t.phi).min(
                            ) for t in list(self.gen_arch.traits.values())])
        else:
            min_fit = 1 - np.atleast_2d(self.gen_arch.traits[
                                                        trt_num].phi).min()

        #then get uneven cmap and cbar-maker (needs to be uneven to
        #give color-resolution to values varying
        #between 1 and the minimum-fitness value assuming all
        #phenotypes are constrained 0<=z<=1, but then also
        #allow gradually deepening reds for fitness values lower
        #than that minimum value), using the min_fit val
        cmap, make_cbar_fn = viz._make_fitness_cmap_and_cbar_maker(
            min_val = min_fit, max_val = 1, cmap = fit_cmap,
                                                trt_num = trt_num)

        #plot the trait phenotype in larger circles first, if trait is not None
        if trt_num is not None:
            #plot the outer (phenotype) circles
            self._plot_phenotype(trait = trt_num, lyr_num = lyr_num,
                land = land, individs = individs, text = False, size = size,
                text_size = text_size, edge_color=edge_color,
                text_color = text_color, cbar = cbar, alpha = alpha,
                zoom_width = zoom_width, x = x, y = y, ticks=ticks,
                mask_rast=mask_rast)
            #make size smaller for the next layer of inner (fitness) circles
            size = round(0.4*size)
            # get sizes for all individuals' inner-circle fitness points, if
            # trt_num is not None
            if trt_num is not None:
                size = size * (1 - ((np.array([*w.values(
                                            )]) - min_fit) / (1 - min_fit)))

        self._plot(lyr_num=lyr_num, land=land, hide_land=True,
                   individs=individs, text=text, color=list(w.values()),
                   pt_cmap=cmap, size=size, edge_color=edge_color,
                   text_color=text_color, cbar=cbar, text_size=text_size,
                   alpha=alpha, zoom_width=zoom_width, x=x, y=y, ticks=ticks,
                   mask_rast=mask_rast)

        #plot phenotype text (works only if plotting a specific trait)
        if phenotype_text and trt_num is not None:
            for ind in self.values():
                plt.text(ind.x-0.5, ind.y-0.5, '%0.2f' % ind.z[trt_num],
                    color = phenotype_text_color, size = text_size)

        #plot fitness text
        if fitness_text:
            offset_from_phenotype_text = 0.001*max(self._land_dim)
            for ind in self.values():
                plt.text(ind.x-0.5+offset_from_phenotype_text,
                         ind.y-0.5+offset_from_phenotype_text,
                         '%0.2f' % ind.fit, color = fitness_text_color,
                         size = text_size)

        #and make a colorbar for the fitness values 
        if fitness_cbar:
            viz._make_fitness_cbar(make_cbar_fn, min_fit)

    #method to plot a species' allele frequencies
    def _plot_allele_frequencies(self):
        if self.gen_arch is None:
            print(("Species._plot_allele_frequencies is not valid for "
                   "species without genomes.\n"))
        else:
            self.gen_arch._plot_allele_frequencies(self)

    # method for plotting a histogram of the current fitness values
    def _plot_hist_fitness(self):
        plt.hist(list(self._calc_fitness()))
        plt.xlabel('Fitness')
        plt.ylabel('Count')

    # method for plotting the movement surface (in various formats)
    def _plot_direction_surface(self, land, surf_type, style, x=None, y=None,
                                zoom_width=None, scale_fact=4.5,
                                color='black', cbar=True, ticks=None,
                                cmap='plasma', mask_rast=None):
        # get the correct surface
        if surf_type == 'move':
            surf = self._move_surf
        elif surf_type == 'disp':
            surf == self._disp_surf

        # get all x's and y's, if x and y are None
        if x is None and y is None:
            x = [*range(land.dim[0])]
            y = [*range(land.dim[1])]
        else:
            x = [x]
            y = [y]
        #check if the surface is none
        if surf is None:
            print(('Function not valid for a Species with no '
                   '_%sSurface.') % ({
                'move': 'Movement', 'disp': 'Dispersal'}[surf_type]))
            return
        elif style not in ['hist', 'chist', 'vect', 'cdraws']:
            print(("The 'style' argument must take one of the "
                   "following values: 'hist', 'chist', "
                   "'vect', 'cdraws'"))
            return
        elif style == 'hist':
            x = x[0]
            y = y[0]
            plt.hist(r.choice(surf.surf[y,x,:], size = 10000, replace = True),
                     bins=100, density=True, alpha=0.5, color=color)

        else:
            #display the movement-surface raster
            lyr_num = surf.lyr_num
            land[lyr_num]._plot(zoom_width = zoom_width, x=np.mean(x),
                                y=np.mean(y), ticks=ticks, cmap=cmap,
                                mask_rast=mask_rast)

            if style == 'chist':
                for x_val in x:
                    for y_val in y:
                        v, a = np.histogram(r.choice(surf.surf[y_val,
                                                               x_val, :],
                                                     replace=True,
                                                     size=7500), bins=15)
                        v = v / float(v.sum())
                        a = [(a[n] + a[n + 1]) / 2 for n in range(len(a) - 1)]
                        xs = [np.cos(a[n]) * 0.75 for n in range(len(a))]
                        ys = [np.sin(a[n]) * 0.75 for n in range(len(a))]
                        xs = np.array(xs) * v * scale_fact
                        ys = np.array(ys) * v * scale_fact
                        [plt.plot((x_val + 0.5, (x_val + 0.5 + xs[n])),
                                  (y_val + 0.5, (y_val + 0.5 + ys[n])),
                                  linewidth=2,
                                  color=color) for n in range(len(xs))]

            elif style == 'cdraws':
                for x_val in x:
                    for y_val in y:
                        pts = [(np.cos(a), np.sin(a)) for a in r.choice(
                                surf.surf[y_val, x_val, :], size=1000,
                                replace=True)]
                        plt.scatter([pt[0] * 0.5 + x_val + 0.5 for pt in pts],
                                    [pt[1] * 0.5 + y_val + 0.5 for pt in pts],
                                    color=color, alpha=0.1, marker='.')

            elif style == 'vect':
                def plot_one_cell(x, y):
                    # draw sample of angles from the Gaussian KDE
                    #representing the von mises mixture distribution (KDE)
                    samp = surf.surf[y,x,:]
                    # create lists of the x and y (i.e. cos and sin)
                    #components of each angle in the sample
                    x_vects = np.cos(samp)
                    y_vects = np.sin(samp)
                    # define the dx and dy distances used to the
                    #position the arrowhead (divide by sqrt(2)/2, to 
                    #scale to the size of half of the diagonal of a cell)
                    dx = np.mean(x_vects) / np.sqrt(2)
                    dy = np.mean(y_vects) / np.sqrt(2)
                    # now plot the arrow
                    plt.arrow(x + 0.5, y + 0.5, dx, dy, alpha=0.75,
                              color=color, head_width=0.24, head_length=0.32)

                # call the internally defined function as a nested list
                #comprehension for all raster cells, which I believe
                #should do its best to vectorize the whole operation
                [[plot_one_cell(j, i) for i in range(
                    surf.surf.shape[0])] for j in range(
                                    surf.surf.shape[1])]


    # method for plotting a species' population pyramid
    def _plot_demographic_pyramid(self):
        #make dict of female and male colors
        col_dict = {-1: 'cyan', 1: 'pink'}
        #create a figure
        fig = plt.figure()
        #variables to grab the max count
        max_count = 0
        #for each sex
        for sex_val in [*col_dict]:
            #get a counter
            counts = C([ind.age for ind in self.values(
                                            ) if ind.sex == int(sex_val < 0)])
            #grab the ages from it
            ages = [*counts]
            #and grab the counts from it (multiplying by -1 for females,
            #to set one sex on either side of x=0, for the pyramid)
            counts = [sex_val*count for count in counts.values()]
            #update the max_count var
            max_count = max(max(counts), max_count)
            #then create the horizontal barplot
            plt.barh(ages, counts, color = col_dict[sex_val])
        #use max_count to set the x-limits and y-limits
        plt.xlim((-1*max_count-2, max_count+2))
        #set the axis labels
        plt.xlabel('Species (individuals)')
        plt.ylabel('Age (timesteps)')
        #then set the xlabels to positive numbers on either side
        locs, labels = plt.xticks()
        plt.xticks(locs, [str(int(loc)) for loc in np.abs(locs)])
        #add sex symbols as title
        plt.suptitle('\u2642%s\u2640' % ''.join([' ' for _ in range(20)]),
                                                                    size = 30)
        #show it
        plt.show()

    def _plot_pop_growth(self):
        T = range(len(self.Nt))
        x0 = self.Nt[0] / self.K.sum()
        plt.plot(T, [_calc_logistic_soln(x0, self.R,
                                t) * self.K.sum() for t in T], color='red')
        plt.plot(T, self.Nt, color='blue')
        plt.xlabel('t')
        plt.ylabel('N(t)')

    def _plot_demographic_changes(self):
        if self._changer is None:
            print(("Species._plot_demographic_changes is not valid "
                "for species with no _SpeciesChanger object.\n"))
        else:
            self._changer._plot_dem_changes(self)

    def _plot_example_recombinant_genome(self):
        assert self.gen_arch is not None, ("This species does not have "
                                            "genomes, so it cannot be used "
                                            "to plot an example recombinant "
                                            "genome.")

        recomb_paths = self.gen_arch._recomb_paths._get_paths(2)
        idxs = (0, 2)
        mock_spp = {}
        for idx in idxs:
            new_genome = np.hstack([(np.ones((self.gen_arch.L,
                                              1)) * n ) + idx for n in (1,2)])
            mock_spp[idx] = Individual(idx=idx, x=0, y=0,
                                       new_genome=new_genome)

        recomb_genome = _do_mating(spp=mock_spp,
                                   mating_pairs=[idxs],
                                   n_offspring=[1],
                                   recomb_paths=recomb_paths)[0][0]
        #recomb_genome = recomb_genome+1
        #recomb_genome[:,1] = recomb_genome[:,1] * 5
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(np.repeat(recomb_genome,
                            int(0.05 * self.gen_arch.L),
                            axis=1),
                  cmap='terrain')
        ax.set_title("gamete 0     gamete1")
        ax.set_xticks([])
        ax.set_ylabel('locus')
        plt.show()
        return(recomb_paths, mock_spp, recomb_genome)


    def _plot_stat(self, stat):
        if self._stats_collector is None:
            print(("Species._plot_stat is not valid "
                "for species with no _StatsCollector object.\n"))

        else:
            self._stats_collector._plot_stat(stat, spp_name = self.name)


        ################
        #public methods#
        ################

    def _write_pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#function to be called when a Species is initiated,
#which uses the params.species[<spp_num>].init.['K_<>'] parameters 
#to make the initial carrying-capacity raster, and also
#sets the Species' K_layer and K_factor attributes
def _make_K(spp, land, K_layer, K_factor):
    #make sure we find only a single layer with the name specified by K_layer
    K_layer = [lyr for lyr in land.values() if lyr.name == K_layer]
    assert len(K_layer) == 1, ("The K_layer parameter should point to"
        "a single Layer, but instead %i Layers were found.") % len(K_layer)
    #grab the identified layer
    K_layer = K_layer[0]
    #set the Species' K_layer and K_factor attributes
    spp.K_layer = K_layer.idx
    spp.K_factor = K_factor
    #add this Species to this Layer's ._is_K attribute
    #(which will be used to update Species.K if this Layer undergoes
    #any kandscape changes)
    K_layer._is_K.append(spp.idx)
    #now calculate and set the K raster
    spp._set_K(land)


def _make_species(land, name, idx, spp_params, burn=False, verbose=False):
    #get spp's intializing params
    init_params = deepcopy(spp_params.init)

    # print verbose output
    if verbose:
        print('\t\tMAKING SPECIES %s...\n' % name)

    #if this species should have genomes, create the genomic architecture
    if 'gen_arch' in spp_params.keys():
        # print verbose output
        if verbose:
            print('\t\t\tmaking genomic architecture...\n')
        g_params = spp_params.gen_arch
        #make genomic_architecture
        gen_arch = _make_genomic_architecture(spp_params = spp_params,
                                                                land = land)
    else:
        gen_arch = None

    # print verbose output
    if verbose:
        print('\t\t\tmaking individuals...\n')
    #make individs
    N = init_params.pop('N')
    #create an ordered dictionary to hold the individuals, and fill it up
    inds = OD()
    for ind_idx in range(N):
        # use individual.create_individual to simulate individuals
        #and add them to the species
        ind = _make_individual(idx=ind_idx, offspring=False,
                               dim=land.dim, genomic_architecture=gen_arch,
                               burn=burn)
        inds[ind_idx] = ind

    #create the species from those individuals
    spp = Species(name = name, idx = idx, inds = inds, land = land,
                     spp_params = spp_params, genomic_architecture=gen_arch)

    #use the remaining init_params to set the carrying-capacity raster (K)
    _make_K(spp, land, **init_params)
    #set initial environment values
    spp._set_e(land)
    #set initial coords and cells
    spp._set_coords_and_cells()
    #set the kd_tree
    spp._set_kd_tree()

    #set phenotypes, if the species has genomes
    #if spp.gen_arch is not None and not burn:
    #    spp._set_z()

    #make density_grid
    spp._set_dens_grids(land)

    #make movement surface, if needed
    if spp._move:
        if 'move_surf' in spp_params.movement.keys():
            if verbose:
                print(('\t\t\tmaking movement surface...\n'
                       '\t\t\t\t[can take a bit]\n'))
            ms_params = deepcopy(spp_params.movement.move_surf)
            #grab the lyr number for the lyr that the 
            #movement surface is to be based on
            move_surf_lyr = ms_params.pop('layer')
            move_surf_lyr_num = [k for k,v in land.items(
                                            ) if v.name == move_surf_lyr]
            assert len(move_surf_lyr_num) == 1, ("Expected to find only a "
                "single Layer with the name provided for the "
                "_ConductanceSurface,"
                " but instead found %i") % len(move_surf_lyr_num)
            move_surf_lyr_num = move_surf_lyr_num[0]
            #make the movement surface and set it as the spp's
            #move_surf attribute
            spp._move_surf= spt._ConductanceSurface(land[move_surf_lyr_num],
                                                                **ms_params)
    #make dispersal surface, if needed
    if 'disp_surf' in spp_params.movement.keys():
        # print verbose output
        if verbose:
            print(('\t\t\tmaking dispersal surface...\n'
                   '\t\t\t\t[can take a bit]\n'))
        ds_params = deepcopy(spp_params.movement.disp_surf)
        #grab the lyr number for the lyr that the 
        #dispersal surface is to be based on
        disp_surf_lyr = ds_params.pop('layer')
        disp_surf_lyr_num = [k for k,v in land.items(
                                        ) if v.name == disp_surf_lyr]
        assert len(disp_surf_lyr_num) == 1, ("Expected to find only a "
            "single Layer with the name provided for the "
            "_ConductanceSurface, "
            "but instead found %i") % len(disp_surf_lyr_num)
        disp_surf_lyr_num = disp_surf_lyr_num[0]
        #make the dispersal surface and set it as the spp's
        #disp_surf attribute
        spp._disp_surf = spt._ConductanceSurface(land[disp_surf_lyr_num],
                                                                **ds_params)

    #if this species has changes parameterized, or if not but it has
    #either a MovementSurf or a DispersalSurf based on a Layer that
    #will undergo landscape change, then create a _SpeciesChanger object for it
    if ('change' in spp_params.keys()
        or (spp._move_surf is not None
        and land._changer is not None
        and spp._move_surf.lyr_num in land._changer.change_info.keys())
        or (spp._disp_surf is not None
        and land._changer is not None
        and spp._disp_surf.lyr_num in land._changer.change_info.keys())):
        # print verbose output
        if verbose:
            print(('\t\t\tsetting up species changes...\n'
                   '\t\t\t\t[can take a while,\n\t\t\t\t if movement or '
                   'dispersal\n\t\t\t\t surfaces will change]\n'))
        #grab the change params (or None, if
        if 'change' in spp_params.keys():
            ch_params = spp_params.change
        else:
            ch_params = None
        #make _SpeciesChanger and set it to the spp's changer attribute
        spp._changer = _SpeciesChanger(spp, ch_params, land = land)

    return spp


# function for reading in a pickled spp
def read_pickled_spp(filename):
    import cPickle
    with open(filename, 'rb') as f:
        spp = cPickle.load(f)
    return spp


