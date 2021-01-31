#!/usr/bin/python
# species.py

# flake8: noqa

'''
Defines the Species class, with its associated methods and supporting functions
'''

#geonomics imports
from geonomics.utils import viz, spatial as spt
from geonomics.structs.genome import (_make_genomic_architecture,
                                      _check_mutation_rates,
                                      _make_starting_mutations,
                                      _get_lineage_dicts,
                                      _get_lineage_dicts_one_tree,
                                      _get_treenums,
                                      _calc_lineage_stat)
from geonomics.structs.landscape import Layer
from geonomics.structs.individual import Individual, _make_individual
from geonomics.ops.movement import _do_movement, _do_dispersal
from geonomics.ops.mating import _find_mates, _draw_n_births, _do_mating
from geonomics.ops.selection import _calc_fitness
from geonomics.ops.mutation import (_do_mutation,
                                    _calc_estimated_total_mutations)
from geonomics.ops.demography import _do_pop_dynamics, _calc_logistic_soln
from geonomics.ops.change import _SpeciesChanger
from geonomics.sim import burnin

#other imports
import numpy as np
from numpy import random as r
import pandas as pd
import random
import matplotlib as mpl
viz._check_display()
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import CCA
from scipy.stats.distributions import norm
from math import pi as math_pi
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
    """
    Representation of the total simulated population of a species.

    Organized as a collections.OrderedDict of serial integer-keyed Individuals.
    Those serial indices continually increment through model time as new
    Individuals (i.e. offspring) as produced, such that no two Individuals
    within the full history of a simulation will ever have the same index
    number.

    Because the Species class inherits from `collections.OrderedDict`,
    Individuals can be indexed out using their index-number keys
    (e.g. Species[<idx>]).

    All Species within a simulation are stored as serial integer-keyed values
    within the Community dict (which in turn is stored at the Model's
    'mod.comm' attribute).

    Attributes
    ----------

        NOTE: For more detail, see the documentation for the parameters that
              correspond to many of the following attributes.

        b:
            The Species' intrinsic birth rate, expressed as the probability (in
            a Bernoulli random draw) that an otherwise eligible potential
            mating pair will successfully reproduce

        burned:
            A bool flag indicating whether or not the Species has been
            burned in.

        coords:
            A 2d numpy array, length(Species) x 2, composed of columns
            containing the current x and y locations of all Individuals.
            Row order is equivalent to the order of individuals in the Species'
            OrderedDict.

        cells:
            A 2d numpy array, length(Species) x 2, composed of columns
            containing the current j and i cell numbers of all Individuals.
            Row order is equivalent to the order of individuals in the Species'
            OrderedDict. (Note the order of the columns! These are not in the
            i,j order that would be necessary in order use them to subject a
            Layer.rast object!)

        choose_nearest_mate:
            A bool flag indicating whether or not Individuals should 
            always choose their nearest neighbors as their potential mates

        d_max:
            The maximum probability of death that can be assigned to any
            Individual

        d_min:
            The minimum probability of death that can be assigned to any
            Individual

        density_grid_window_width:
            The width of the sliding window that is used by the Species'
            _DensityGridStack object when estimating the Species'
            population-density array.

        direction_distr_[mu/kappa]:
            The $\mu$ and $\kappa$ parameters of the von Mises distribution
            from which the direction components of Individuals' movement
            vectors will be drawn. (These parameters are only utilized if the
            Species is parameterized without a _MovementSurface.)

        dispersal_distance_distr:
            A string indicating which distribution to use as
            the Species' dispersal-distance distribution ('lognormal',
            'wald', or 'levy')

        dispersal_distance_distr_param[1/2]:
            The values of the first and second parameters of the Species'
            dispersal-distance distribution (named abstractly like this
            because they could parameterize a Lognormal, Wald,
            or Levy distribution, depending on the value of
            Species.dispersal_distance_distr)

        extinct:
            A bool flag indicating whether or not the Species has gone extinct.

        idx:
            Index number of the Species (i.e. its key within the Community dict)

        inverse_dist_mating:
            A bool flag indicating whether or not the probabilities that an
            Individual chooses any of its neighbors (i.e. other Individuals
            within `Species.mating_radius` distance of them) as their
            potential mates should vary inversely with the distances between
            the Individual and each of those neighbors. If False, any neighbor
            could be chosen with equal probability. 

        gen_arch:
            The Species' genomic architecture (as a GenomicArchitecture object)

        K:
            The Species' current carrying-capacity raster. This is a numpy array
            of shape `Landscape.dim`, with each value 
            indicating the carrying capacity of that cell in the Landscape.
            It is generated as the product of the Species' K_factor
            and the raster array of the Layer provided as the basis
            of the Species' carrying capacity
            (i.e. `spp.K_factor * mod.land[spp.K_layer].rast`).
            This array will change during a simulation if the corresponding
            Layer is parameterized to undergo environmental change.

        K_factor:
            The factor by which to multiply the Species' chosen
            carrying capacity Layer (mod.land[spp.K_layer]) in order to
            generate the Species' carrying-capacity raster (spp.K).

        K_layer:
            The index number of the Layer whose raster will be multiplied by
            spp.K_factor in order to generate the Species' carrying-capacity
            raster (spp.K).

        mating_radius:
            The mating radius, expressed in Landscape cell widths, from within
            which an Individual can choose a potential mate. This will only be
            used if `Species.choose_nearest_mate` is False.

        max_age:
            The maximum age, in time steps, that any Individual of the Species
            may reach (after which point they will be culled from the
            population by imposed mortality)

        max_ind_idx:
            The maximum Individual index number that has currently been used

        move:
            A bool flag indicating whether or not Individuals of the Species
            can move after their initial dispersal as offspring.

        movement_distance_distr:
            A string indicating which distribution to use as
            the Species' movement-distance distribution ('lognormal',
            'wald', or 'levy')

        movement_distance_distr_param[1/2]:
            The values of the first and second parameters of the Species'
            movement-distance distribution (named abstractly like this
            because they could parameterize a Lognormal, Wald,
            or Levy distribution, depending on the value of
            Species.movement_distance_distr)

        mutate:
            A bool flag indicating whether or not the Species' simulation
            involves mutation.

        mut_log:
            A bool flag indicating whether or not a mutation log should be
            written for the Species.

        N:
            The Species' current population-density raster.
            This is a numpy array of shape `Landscape.dim`, with each value 
            indicating the population density of that cell in the Landscape,
            as estimated by the Species' _DensityGridStack instance. 

        Nt:
            A chronologically ordered list, starting from the first timestep of
            the burn-in, containing the Species' total population size at each
            time step.

        n_births:
            A chronologically ordered list, starting from the first timestep of
            the burn-in, containing the Species' number of births at each
            time step.

        n_births_distr_lambda:
            If `Species.n_births_fixed` is False, then this value serves as the 
            $\lambda$ parameter of the Poisson distribution from which is drawn
            the number of offspring a mating pair will produce.
            If `Species.n_births_fixed` is True, then this is the fixed number
            of offspring that each successful mating pair will produce.

        n_births_fixed:
            A bool flag indicating whether or not the number of births should
            be fixed at `Species.n_births_distr_lambda`

        n_deaths:
            A chronologically ordered list, starting from the first timestep of
            the burn-in, containing the Species' number of deaths at each
            time step.

        name:
            The string name of the Species

        R:
            The intrinsic growth rate of the Species, to be used as 'R' in the
            stochastic, spatialized logistic growth model that controls
            population dynamics

        repro_age:
            The age, in time steps, of first reproduction

        selection:
            A bool flag indicating whether or not the Species' simulation
            involves selection.

        sex:
            A bool flag indicating whether or not this Species is sexed
            (i.e. whether or not the eligibility of potential mating pairs
            should be basd on whether or not the Individuals in the potential
            pair are of opposite sexes)

        sex_ratio:
            The sex ratio of the Species. This value is expressed as
            the proportion of all offspring that are males, such that it can be
            easily used as the probability parameter for the Bernoulli draw
            of an offspring's sex. (Note that this value is derived from
            the 'sex_ratio' parameter provided in the Model's parameters file, 
            but unlike in the parameters file this value is
            not expressed in as a ratio of males to females.)

        start_N:
            The number of Individuals to be generated for the starting
            population of the Species (i.e. the population size at the
            beginning of the burn-in)

        t:
            The index of the last completed time step run for this Species.
            This attribute holds the value -1 until the Species has been burned
            in, at which point it begins incrementing, Pythonically, from 0
            (such that a Species with t == 999 has been run for 1000 time
            steps).

    """
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
            self._tc_sorted_and_simplified = False
        else:
            self._tc = None
            self._tc_sorted_and_simplified = None

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

        # create the burn-in spatial counter
        self._burnin_spat_tester = burnin.SpatialTester(self)


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

    # method for running the spatial burn-in test
    def _do_spatial_burnin_test(self, num_timesteps_back):
        # first update the counts and diffs arrays
        # (which count the number of individuals within each landscape cell,
        # and the diff in that number between last check and this one)
        self._burnin_spat_tester.update(self)
        result = self._burnin_spat_tester.run_test(num_timesteps_back)
        return result

    #method to move all individuals simultaneously, and sample
    #their new locations' environment
    def _do_movement(self, land):
        _do_movement(self)
        self._set_e(land)
        self._set_coords_and_cells()

    #function for finding all the mating pairs in a species
    def _find_mating_pairs(self):
        mating_pairs = _find_mates(self, sex=self.sex,
                                   choose_nearest=self.choose_nearest_mate,
                                   inverse_dist_mating=self.inverse_dist_mating)
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
        offspring_keys = list(range(next_offspring_key,
                                    next_offspring_key + total_births))[::-1]
        #update self.max_ind_idx
        if len(offspring_keys) > 0:
            self.max_ind_idx = offspring_keys[0]

        #copy the keys, for use in mutation.do_mutation()
        keys_list = [*offspring_keys]

        if not burn and self.gen_arch is not None:
            recomb_keys = [*r.randint(low=0,
                                      high=self.gen_arch.recombinations._n,
                                      size=total_births*2)]
            # NOTE: this gives us a hierarchically organized list,
            # containing both the offspring's new genomes and their
            # segment information (to be added to the
            # tskit.TableCollection.edges table), organized as a sublist for
            # each pair, then a subtuple within that, containing a new genome
            # and then a 2-tuple of segments info for each of the new
            # offspring's homologues
            genomes_and_segs = _do_mating(self, mating_pairs,
                                          n_births, recomb_keys)

        for n_pair, pair in enumerate(mating_pairs):

            parent_midpoint_x = (self[pair[0]].x + self[pair[1]].x)/2
            parent_midpoint_y = (self[pair[0]].y + self[pair[1]].y)/2

            n_offspring = n_births[n_pair]

            for n in range(n_offspring):

                #get the next offspring_key
                offspring_key = offspring_keys.pop()

                offspring_x, offspring_y = _do_dispersal(
                    self, parent_midpoint_x, parent_midpoint_y,
                    self.dispersal_distance_distr_param1,
                    self.dispersal_distance_distr_param2)

                #set the age to 0
                age = 0

                #set the sex correctly
                if self.sex:
                    sex = r.binomial(1, self.sex_ratio)
                else:
                    sex = None

                #set the new_genome correctly
                if self.gen_arch is not None:
                    if not burn:
                        new_genome = genomes_and_segs[n_pair][n][0]
                        genome_segs = genomes_and_segs[n_pair][n][1]
                    else:
                        new_genome=None
                else:
                    new_genome = None

                #create the new individual
                self[offspring_key] = Individual(idx=offspring_key, age=age,
                                                 new_genome=new_genome,
                                                 x=offspring_x, y=offspring_y,
                                                 sex=sex)

                #set new individual's phenotype (won't be set
                #during burn-in, because no genomes assigned;
                #won't be set if the species has no gen_arch)
                if (self.gen_arch is not None
                    and self.gen_arch.traits is not None
                    and not burn):
                    self[offspring_key]._set_z(self.gen_arch)

                # during the main phase, for species with genomes,
                # update the tskit tables
                if (self.gen_arch is not None
                    and not burn):
                    loc = [offspring_x, offspring_y]
                    if self.gen_arch.traits is not None:
                        loc = (loc + self[offspring_key].z +
                            [self[offspring_key].fit])
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
                    # the 'flags' column vals to 1
                    # (to indicate they're real individs, not msprime-derived)
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
                        # TODO: DOES THIS LINE THROW AN ERROR B/C IT LOOPS
                        # OVER >2 seg_sets, SO INDEXING FOR homol = >1 IS
                        # A KEY ERROR FOR A DICT WITH KEYS 0 AND 1?
                        child=self[offspring_key]._nodes_tab_ids[homol]
                        ) for homol, seg_set in enumerate(
                        genome_segs) for seg in [*seg_set]]


        # sample all individuals' environment values, to initiate for offspring
        self._set_e(land)
        self._set_coords_and_cells()

        # do mutation if necessary
        if self.mutate and not burn:
             _do_mutation(keys_list, self, log = self.mut_log)


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

    # method to set the species' genomes to all zeros (after burn-in),
    # if the species has any nonneutral loci or has non-zero nonneutral
    # mutation rates
    def _set_null_genomes(self):
        if (len(self.gen_arch.nonneut_loci) > 0 or
            (self.gen_arch.traits is not None and
             len([trt.mu > 0 for trt in self.gen_arch.traits]) > 0) or
            self.gen_arch.mu_delet > 0):
            [ind._set_g(np.zeros(
                (len(self.gen_arch.nonneut_loci),
                    self.gen_arch.x))) for ind in self.values()]

    # add new row to the individuals' numpy arrays, for a given mutation locus
    def _add_new_locus(self, idx, locus):
        # insert a row of ones in each individual's genotype array 
        [ind._add_new_locus(idx) for ind in self.values()]

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
    def _set_z(self):
        [ind._set_z(self.gen_arch) for ind in self.values()];

    #method for setting an individual's phenotype
    def _set_z_individ(self, individ):
        self[individ]._set_z(self.gen_arch)

    #method to set the individuals' fitness attributes
    def _set_fit(self, fit):
        [ind._set_fit(f) for ind, f in zip(self.values(), fit)];

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


    # method to set the species' genomes, check for adequate mutable loci,
    # and fill the tskit.TableCollection's tables
    # (to be called after the model has burned in)
    def _set_genomes_and_tables(self, burn_T, T):
        # set the species' neutral and non-netural loci,
        # and set genomes to all zeros
        self._set_null_genomes()

        # use mean n_births at tail end of burn-in to estimate number of
        # mutations for this species, then warn the user if there is or may be
        # inadequate space for the parameterized mutation rate (because only
        # using infinite-sites mutation)
        est_tot_muts = _calc_estimated_total_mutations(self, burn_T, T)

        # check whether there are adequate mutable loci for this species
        _check_mutation_rates(self.gen_arch, est_tot_muts, burn_T, T)

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
    def _sort_simplify_table_collection(self, check_nodes=False,
                                        check_individuals=False):
        # sort the TableCollection
        self._tc.sort()
        # get an array of all the current individuals' nodes,
        # the nodes for which the tables will be simplified
        curr_nodes = self._get_nodes()

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
        if check_individuals:
            print('b4', len(meta_b4), 'af', len(ind_meta))
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
                print('...\ngnx id:', ind, '\ttable id:', idx)
                print('gnx node ids:  ', [*self[ind]._nodes_tab_ids.values()])
                print('table node ids:', [*subdf.index])
                node_id_vals = [*subdf.index]
                row_idx = np.where(new_ids_from_tables[:,0] == ind)[0][0]
                new_ids_from_tables[row_idx, 1:] = node_id_vals
            # now make sure tabls are identical
            assert np.all(new_ids_from_tables == new_ids), ("All the node "
                                                            "ids don't match!")
            print('ALL NODE IDS APPEAR CORRECTLY REASSIGNED\n')

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
            print(z0_af)
            z0_check = [z0 == z0_b4[idx] for idx, z0 in z0_af.items()]
            for idx in [*self]:
                print(idx, ' b4 ', z0_b4[idx])
                print(idx, ' af ', z0_af[idx])
                print('..................')
            assert np.all(z0_check), ('phenotypes are not the same! '
                                      ' (MAKE SURE YOU USED A MODEL WITH '
                                      ' AT LEAST ONE TRAIT!)')
            print('ALL INDIVIDUAL IDS APPEAR CORRECTLY REASSIGNED\n')

        # set the sorted_and_simplified flag to True
        self._tc_sorted_and_simplified = True


    # get the nodes-table IDs for all individuals
    def _get_nodes(self, individs=None):
        if individs is None:
            individs = [*self]
        nodes = np.hstack([[*self[i]._nodes_tab_ids.values(
                                                        )] for i in individs])
        return nodes


    ##########################################
    # FUNCTIONS FROM TRACK_SPATIAL_PEDIGREE.PY
    ##########################################

    # for a sample of nodes and a sample of loci, get the
    # loci-sequentially ordered dict of lineage dicts
    # (containing each node in a child node's lineage as the
    # key and its birth-time and birth-location as a length-2 list of values
    def _get_lineage_dicts(self, loci, nodes=None, drop_before_sim=True,
                           time_before_present=True,
                           use_individs_curr_pos=True,
                           max_time_ago=None,
                           min_time_ago=None):
        if nodes is None:
            #sort and simplify the TableCollection, if needed
            if not self._tc_sorted_and_simplified:
                self._sort_simplify_table_collection()
            nodes = self._get_nodes()
        lin_dicts = _get_lineage_dicts(self, nodes, loci, t_curr=self.t,
                                       drop_before_sim=drop_before_sim,
                                       time_before_present=time_before_present,
                                       max_time_ago=max_time_ago,
                                       min_time_ago=min_time_ago)
        # if requested, put the current individuals' current 
        # positions in the lineage dicts, rather than using the birth
        # positions (which are in there by default because those
        # are the positions that are stored within the tskit tables)
        if use_individs_curr_pos:
            for locus in lin_dicts.keys():
                for curr_node_id in lin_dicts[locus].keys():
                    if len(lin_dicts[locus][curr_node_id]) > 0:
                        tc_individs_id = self._tc.nodes[curr_node_id].individual
                        individ = [ind for ind in self.values() if
                                   ind._individuals_tab_id == tc_individs_id]
                        assert len(individ) == 1, ("Found multiple individuals "
                                                   "with the same individuals id "
                                                   "in the tskit tables.")
                        individ = individ[0]
                        curr_pos = np.array([individ.x, individ.y])
                        birth_t = lin_dicts[locus][curr_node_id][curr_node_id][0]
                        lin_dicts[locus][curr_node_id][curr_node_id] = (birth_t,
                                                                        curr_pos)
        return lin_dicts


    # check whether specified individuals have coalesced at specified loci
    # (within the forward-time Geonomics run, of course, since by design
    # everything coalesces in the ms-generated pre-simulation pedigree)
    def _check_coalescence(self, individs=None, loci=None, all_loci=False):
        #sort and simplify the TableCollection, if needed
        if not self._tc_sorted_and_simplified:
            self._sort_simplify_table_collection()
        if individs is None:
            individs = [*self]
        # get the node IDs
        nodes = self._get_nodes(individs=individs)
        if loci is None:
            loci = [*range(self.gen_arch.L)]
        # NOTE: if running into memory issues when getting the lineage dicts
        # all at once for simulations with large pop sizes and large nums of
        # loci then should try instead getting them 1 by 1 in the loop below
        lin_dicts = self._get_lineage_dicts(loci, nodes=nodes)
        result = {}
        for loc, lin_dict in lin_dicts.items():
            # get the number of unique oldest nodes across the lineage dict
            oldest_nodes = [[*lin_dict[node].keys(
                                            )][-1] for node in lin_dict.keys()]
            n_unique_nodes = len(np.unique(oldest_nodes))
            result[loc] = n_unique_nodes == 1
        # check coalescence across all loci, if requested
        if all_loci:
            result = np.all([*result.values()])
        return result


    def _calc_lineage_stats(self, individs=None, nodes=None, loci=None,
                            stats=['dir', 'dist', 'time', 'speed'],
                            use_individs_curr_pos=True,
                            max_time_ago=None, min_time_ago=None):
        '''
        Calculate stats for the lineages of a given set of nodes and loci;
        returns dict of struct: {k=stat, v={k=loc, v=[val_node1 ... val_node_N]}}
        '''
        # get all nodes for the provided individuals, or for all individuals,
        # if nodes IDs not provided
        if nodes is None:
            #sort and simplify the TableCollection, if needed
            if not self._tc_sorted_and_simplified:
                self._sort_simplify_table_collection()
            nodes = self._get_nodes(individs)
        # get all loci, if loci not provided
        if loci is None:
            loci = [*range(self.gen_arch.L)]
        lin_dicts = self._get_lineage_dicts(loci, nodes=nodes,
                use_individs_curr_pos=use_individs_curr_pos,
                                            max_time_ago=max_time_ago,
                                            min_time_ago=min_time_ago)
        stats = {stat: {} for stat in stats}
        for stat in stats:
            for loc in loci:
                loc_stat_list = []
                for node in nodes:
                    loc_stat_list.append(_calc_lineage_stat(
                                            lin_dicts[loc][node], stat=stat))
                stats[stat][loc] = loc_stat_list
        return stats


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


    def _get_genotypes(self, loci=None, individs=None, biallelic=True,
                       as_dict=False):
        # sort and simplify the table collection
        # (thus dropping any unnecessary individuals' data in there and also
        # making the tables' structure simler and more predictable)
        self._sort_simplify_table_collection()
        # then get the TreeSequence
        ts = self._tc.tree_sequence()

        # make sure as_dict and biallelic are True or False
        assert as_dict in [True, False], ("The 'as_dict' argument "
                                          "must be either "
                                          "True or False.")

        assert biallelic in [True, False], ("The 'biallelic' argument "
                                          "must be either "
                                          "True or False.")

        # get the list of loci
        if loci is not None:
            assert np.iterable(loci), ("The 'loci' argument needs "
                                       "either None or "
                                       "an iterable of locus indices.")
        # get the list of individuals
        if individs is None:
            individs = [*self]
        else:
            assert np.iterable(individs), ("The 'individs' argument needs "
                                       "an iterable of individual IDs.")
        # sort the individuals
        # NOTE: makes output data strucutes nicer (since individuals' ids are
        # always in numerical order), and also makes it simpler to ensure that
        # there are no mismatching ids and genotypes in downstream data
        # structures
        individs = np.sort(individs)

        # get the list of the individuals' nodes
        samples_to_keep = np.int64(np.hstack([[*self[ind]._nodes_tab_ids.values(
                                                    )] for ind in individs]))
        assert len(samples_to_keep) == self.gen_arch.x * len(individs), ('Num'
                        'ber of nodes does not match number of individs!')

        # get haplotypes for all samples
        # TODO: np.where(self._tc.nodes.flags)[1] gets only the nodes with
        #       flags of 1, i.e. only 'real' individuals who were not part of
        #       the msprime-sourced pre-simulation fake genealogy
        # TODO: haplotypes returned in order of samples???
        #haps = [np.int8([*hap]) for n, hap in zip(
        #                                    np.where(self._tc.nodes.flags)[0],
        #                                    ts.haplotypes()) if n in samples]
        haps_dict = {s:h for s, h in zip(ts.samples(), ts.haplotypes())}
        haps = [np.int8([*haps_dict[s]]) for s in samples_to_keep]

        # get the genotypes by combining each consecutive group
        # of x haplotypes, where x is the ploidy, then vstacking them
        grouped_haps = zip(*[haps[i::self.gen_arch.x] for i in range(
                                                        self.gen_arch.x)])
        gts = [np.vstack(h).T for h in grouped_haps]

        # subset loci, if needed
        if loci is not None:
            gts = [gt[loci, :] for gt in gts]

        # get mean genotype for each individual, if necessary
        if not biallelic:
            gts = [np.mean(gt, axis=1) for gt in gts]

        # stack into the speciome, if dict not requested
        # (dims are N x L x X, where N=num individs, L=num loci,
        #  and X=ploidy if biallelic=True, else X=1)
        if not as_dict:
            gts = np.stack(gts)

        # or else cast as dict
        else:
            gts = {individ: gt for individ, gt in zip(individs, gts)}

        return gts


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
    def _get_z(self, trait_num=None, individs=None):
        zs = self._get_scalar_attr('z', individs=individs)
        if trait_num is not None:
            zs = zs[:,trait_num]
        return zs

    #convenience method for getting whole species' fitnesses
    def _get_fit(self, individs = None):
        fits = self._get_scalar_attr('fit', individs=individs)
        return fits

    def _calc_fitness(self, trait_num = None, set_fit = True):
        fit = _calc_fitness(self, trait_num = trait_num)
        #set individuals' fitness attributes, if indicated
        if set_fit:
            self._set_fit(fit)
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

    #use the kd_tree attribute to find mating pairs either
    #within the species, if within == True, or between the species
    #and the points provided, if within == False and points is not None
    def _get_mating_pairs(self, within=True, coords=None,
                          choose_nearest=False, inverse_dist_mating=False):
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
        #the spp._get_mating_pairs() method is called!
        self._set_kd_tree()

        # if mating_radius is None, then just use Wright-Fisher
        # style panmixia (draw with replacement a sample of size = Nt*b, where Nt is the 
        # current pop size and b is the birth rate (i.e. mating probability);
        # that sample represents all mating individuals, and each of those
        # individuals randomly chooses its mate
        if self.mating_radius is None:
            # draw a number of mating individuals as a binomial rv
            # with num trials equal to pop size and probability equal to
            # the species' birth rate
            if self.b < 1:
                n_mates = np.random.binomial(n=len(self), p=self.b, size=1)
            else:
                n_mates = len(self)
            tree_inds = self._kd_tree.tree.indices
            # draw 2*n_mates mating individuals, with replacement (as in WF
            # model), then fold into an n_mates x 2 mate-pairs array
            pairs = np.random.choice(tree_inds, replace=True,
                                     size=n_mates*2).reshape((n_mates, 2))
            # get rid of selfing pairs
            pairs = np.array([list(pair) for pair in list(
                                    map(set, pairs)) if len(pair) == 2])
            pairs = np.array(pairs)

        # otherwise, choose mates using mating radius
        else:
            #if neighbors are to be found within the species,
            #set coords to self.coords (otherwise, the coords to
            #find nearest neighbors with should have been provided)
            if within:
                coords = self.coords

            #query the tree to get mating pairs
            pairs = self._kd_tree._get_mating_pairs(coords=coords,
                                                    dist=self.mating_radius,
                                                    choose_nearest=choose_nearest,
                                            inverse_dist_mating=inverse_dist_mating)

            # use the species' birth rate to decide (as bernoulli draws)
            # whether each pair can mate
            can_mate = np.bool8(np.random.binomial(n=1, p=self.b,
                                                   size=pairs.shape[0]))
            pairs = pairs[can_mate, :]
        return pairs


    def _make_gea_df(self, lyr_num=1):
        """

        Makes a dataframe for GEA analysis containing genotypes, coordinates,
        and env values for all individuals.

        NOTE: Currently only accepts one Landscape Layer (i.e. environmental
              variable)

        Parameters
        ----------
        lyr_num : int
            The number of the Landscape Layer from which to extract
            environmental values. Defaults to 1.

        Returns
        -------
        out : pandas.DataFrame
            A DataFrame in which columns are loci and rows are individuals
        """
        #get array of [0|0.5|1] genotypes
        """
        gen = self._get_genotypes()

        #loop to convert binary allele genotypes (0|1) into single digit genotypes (0, 0.5, 1)
        geno_ind = []
        for ind in range(len(gen)):
            geno_loc = []
            for loci in range(len(gen[ind])):
                genotype = gen[ind][loci].mean() # coded as 0/0.5/1
                geno_loc.append(genotype)
            geno_ind.append(geno_loc)
        """
        gts = self._get_genotypes(biallelic=False)

        #convert to dataframe
        gea_df = pd.DataFrame(gts)

        #get environmental data
        env = self._get_e()
        gea_df['env'] = list(env[:,lyr_num])

        #get coords
        coord = self._get_coords()
        gea_df['lat'] = list(coord[:,0])
        gea_df['long'] = list(coord[:,1])

        return gea_df


    def _run_cca(self, trt_num=0, scale=True, plot=True, plot_sd=True,
                 sd=3):
        """
        Runs a canonical correlation analysis (CCA) on the current genetic and
        environmental data for the Species' current population, using the
        indicated Trait. Plots the results and returns them in a dict.

        CCA model formula: genotype ~ env + lat + long

        NOTE: currently only possible to do 3 components/3 variables;
              will add greater functionality in the future

        Parameters
        ----------
        trt_num : int
            The number of the Trait to run the GEA on. Defaults to 0.
        scale : bool
            If True, scales the variable, individual, and locus loadings from -1 to 1 to make them easier
            to visualize. Defaults to True.
        plot : bool
            Whether or not to plot the model. Defaults to True.
        plot_sd : bool
            If True, a standard deviation ellipse is plotted, using the number
            of standard deviations indicated by the argument 'sd'. Defaults to
            True.
        sd : {int, float}
            Number of standard deviations to use for plotting the standard
            deviation ellipse. Defaults to 3.

        Returns
        -------
        out : dict
            A dict of the following key-value pairs:
                "gea_df": a DataFrame of individuals' genotypes, environmental
                          values, and coordinates
                "ind_CCAdf": individual loadings (columns = axes,
                             rows = individuals)
                "loci_CCAdf": locus loadings (columns = axes, rows = loci)
                "var_CCAdf": variable loadings (columns = axes,
                             rows = variable (env, lat, long))
                "trait_loci": list of indexes for loci that underlie the trait
                              for which the CCA was run
        """
        #make DF
        gea_df =self._make_gea_df(lyr_num=self.gen_arch.traits[trt_num].lyr_num)

        #get adaptive loci (used later in plotting)
        trait_loci = self.gen_arch.traits[trt_num].loci #gets trait loci indexes

        #get number of loci
        L = self.gen_arch.L
        #define x and y:
        # pull out genotypes
        Y = gea_df.iloc[:,range(L)]
        # pull out predictor variables
        # TODO: make it possible to choose more vars in future
        X = gea_df[['env','lat','long']]

        #define n_components
        n_components = 3
        #create the CCA model and fit it to the data
        cca = CCA(n_components = n_components)
        cca.fit(Y, X)

        #transform data
        Y_c = cca.transform(Y)

        #create df for individuals
        ind_df = pd.DataFrame(Y_c)
        #naming CCA columns starting at 1
        ind_df.columns = [str(r) for r in range(1, n_components + 1)]

        #create df for loci
        loci_df = pd.DataFrame(cca.x_loadings_)
        loci_df.columns = ind_df.columns #same column names

        #create df for variables
        var_df = pd.DataFrame(cca.y_loadings_)
        var_df.columns = ind_df.columns #same column names

        #make dictionary of dataframes to return
        cca_dict = {'gea_df':gea_df,
                    'ind_df':ind_df,
                    'loci_df':loci_df,
                    'var_df':var_df,
                    'trait_loci':trait_loci}

        # NOTE: right now this plotting loop only works for a maximum of 3 axes;
        #       in the future as more env vars are added, could add more axes
        if plot:
            #Scale dataframes from -1 to 1 if scale = True
            if scale == True:
                rmin = -1
                rmax = 1
                loci_df = rmin + (loci_df - loci_df.values.min()) * (rmax - (rmin)) / (loci_df.values.max() - loci_df.values.min())
                var_df = rmin + (var_df - var_df.values.min()) * (rmax - (rmin)) / (var_df.values.max() - var_df.values.min())
                ind_df = rmin + (ind_df - ind_df.values.min()) * (rmax - (rmin)) / (ind_df.values.max() - ind_df.values.min())


            #get max and mins to set axis later on
            maxdf = max(ind_df.values.max(), var_df.values.max(), loci_df.values.max())
            mindf = abs(min(ind_df.values.min(), var_df.values.min(), loci_df.values.min()))
            axmax = max(maxdf, mindf) * 1.20 #adding a 20% buffer
            
            #set up figure
            fig = plt.figure(figsize=(20, 15))
            for n, cc_axes_pair in enumerate([[1, 2], [1, 3], [2, 3]]):
                #define components for axis
                cc_axis1 = cc_axes_pair[0] #x-axis CC
                cc_axis2 = cc_axes_pair[1] #y-axis CC

                ax = fig.add_subplot(131+n)

                #add center lines
                ax.axhline(y=0, color='lightgray', linestyle='dotted')
                ax.axvline(x=0, color='lightgray', linestyle='dotted')
                
                # set axes range
                plt.xlim(-axmax, axmax)
                plt.ylim(-axmax, axmax)

                #plot neutral SNPs
                ax.scatter(loci_df[str(cc_axis1)], loci_df[str(cc_axis2)],
                           marker = '+', color = 'gray',
                           label = 'Neutral SNPs', alpha = 0.7)
                #plot individuals
                cmap_scat = ax.scatter(ind_df[str(cc_axis1)],
                                       ind_df[str(cc_axis2)],
                                       c = gea_df.env, cmap='viridis',
                                       label = "Individuals", alpha = 0.5)
                #plot adaptive SNPs
                ax.scatter(loci_df[str(cc_axis1)][trait_loci],
                           loci_df[str(cc_axis2)][trait_loci],
                           marker="+", color = "red", s=100,
                           label = 'Adaptive SNPs')

                #plot variable vectors as arrows
                for i in range(var_df.shape[0]):
                    x = var_df[str(cc_axis1)][i]
                    y = var_df[str(cc_axis2)][i]
                    plt.arrow(0, 0, x, y, width = 0.01,
                              head_width = 0.05, color = 'black')
                    sx = 0.1
                    sy = 0.1
                    # this mess is just to arrange the text
                    # next to arrows but it is a WIP
                    if (x < 0 and y < 0):
                       plt.text(x - sx, y - sy, X.columns[i])
                    if (x < 0 and y > 0):
                       plt.text(x - sx, y + sy, X.columns[i])
                    if (x > 0 and y > 0):
                       plt.text(x + sx, y + sy, X.columns[i])
                    if (x > 0 and y < 0):
                       plt.text(x + sx, y - sy, X.columns[i])

                #plot SD ellipse for SNP data if plot_sd = True
                if plot_sd:
                    r1 = np.std(loci_df[str(cc_axis1)]) * sd
                    r2 = np.std(loci_df[str(cc_axis2)]) * sd
                    t = np.linspace(0, 2*math_pi, 100)
                    ax.plot(r1*np.cos(t) , r2*np.sin(t), linestyle = 'dashed',
                            c = 'red', label = str(sd) + " StdDev")

                #EVERYTHING BELOW THIS LINE IN THE LOOP IS JUST PLOT FORMATTING 
                #xlabel is the label at the bottom (vertical axis label) and ylabel is on the left side (horizontal axis label)
                label = ax.set_xlabel('CCA' + str(cc_axis2), fontsize = 9)
                ax.xaxis.set_label_coords(0.5, -0.02)
                label = ax.set_ylabel('CCA' + str(cc_axis1), fontsize = 9)
                ax.yaxis.set_label_coords(-0.02, 0.55)

                #make into box
                ax.set_aspect('equal', adjustable='box')
                
                #Move left y-axis and bottom x-axis to center, passing through (0,0)
                ax.spines['left'].set_position('center')
                ax.spines['bottom'].set_position('center')

                #change colors/remove default axis
                axcolor = 'none'
                ax.spines['left'].set_color(axcolor)
                ax.spines['bottom'].set_color(axcolor)
                ax.spines['right'].set_color(axcolor)
                ax.spines['top'].set_color(axcolor)
                ax.tick_params(axis = 'x', colors=axcolor)
                ax.tick_params(axis = 'y', colors=axcolor)

                ax.yaxis.label.set_color('gray')
                ax.xaxis.label.set_color('gray')

                #legends
                if n == 2:
                    ax.legend(loc = "upper left", bbox_to_anchor=(1.3, 1.05))
                    plt.colorbar(cmap_scat, label="env")
            fig.show()

        return cca_dict


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
    def _plot_genotype(self, locus, lyr_num=None, land=None, individs=None,
                       text=False, size=25, text_size = 9, edge_color='black',
                       text_color='black', cbar=True, alpha=1,
                       by_dominance=False, zoom_width=None, x=None, y=None,
                       ticks=None, mask_rast=None):

        genotypes = self._get_genotypes(loci=[locus], individs=individs,
                                        biallelic=False, as_dict=True)
        if by_dominance == True:
            if locus in self.gen_arch.nonneut_loci:
                dom = self.gen_arch.dom[locus]
                # NOTE: would be simpler to just use np.ceil here...
                #       don't see why I didn't...?
                genotypes = {i:np.clip(gt * (1 + dom),
                                       a_min=None,
                                       a_max=1) for i, gt in genotypes.items()}
            else:
                print(('\n\tWARNING: The by_dominance argument is True, '
                       'but a neutral locus was chosen, '
                       'so the argument was not used.\n'))

        # just assign black, gray, and white (since there's no reason
        # necessarily that these should be mapped to a certain layer, the way
        # phenotype should
        colors = ['#000000', '#808080', '#FFFFFF']

        for n, genotype in enumerate([0.0, 0.5, 1.0]):
            genotype_individs = [i for i, g in genotypes.items(
                                        ) if np.atleast_1d(g)[0] == genotype]
            # will hide the land if this is not the first plot made
            hide_land = n > 0
            # plot if there are any individuals of this genotype
            if len(genotype_individs) >= 1:
                self._plot(lyr_num=lyr_num, land=land, hide_land=hide_land,
                           individs=genotype_individs,
                           text=text, color=colors[n], edge_color=edge_color,
                           text_color=text_color, cbar=cbar,
                           size=size, text_size=text_size, alpha=alpha,
                           zoom_width=zoom_width, x=x, y=y, vmin=0, vmax=1,
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

        z = OD(zip([*self], self._get_z()[:, trait]))
        if individs is not None:
            z = {i: v for i, v in z.items() if i in individs}

        # get the correct cmap for this trait's layer
        pt_cmap = viz._choose_cmap(self.gen_arch.traits[trait].lyr_num)

        points = self._plot(lyr_num = lyr_num, land = land,
                            individs = individs, text = text,
                            color = list(z.values()), pt_cmap = pt_cmap,
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


    # plot the lineage for a given node and locus
    def _plot_gene_flow(self, locus, style, land, nodes=None, individs=None,
                       color=None, phenotype=None, lyr_num=0, jitter=True,
                       alpha=0.5, size=25, add_roots=False):
        assert style in ['lineage', 'vector'], ("The style argument must be "
                                                "given either 'lineage' or "
                                                "'vector' as a value.")
        if nodes is None:
            # get a random selection of n individuals, if individs' value is int
            if isinstance(individs, int):
                individs = np.random.choice([*self], individs, replace=False)
            #sort and simplify the TableCollection, if needed 
            if not self._tc_sorted_and_simplified:
                self._sort_simplify_table_collection()
            nodes = self._get_nodes(individs=individs)

        # sort and simplify the TableCollection if needed
        if not self._tc_sorted_and_simplified:
            self._sort_simplify_table_collection()
        # grab the TableCollection and its TreeSequence
        tc = self._tc
        try:
            ts = tc.tree_sequence()
        except Exception:
            raise Exception(("The species' TableCollection must be sorted and "
                             "simplified before this method is called."))

        node_curr_locs = [[i.x, i.y] for n in nodes for i in self.values(
                        ) if n in i._nodes_tab_ids.values()]

        # get the tree for this locus
        tree = ts.aslist()[_get_treenums(ts, [locus])[0]]
        # get the lineage_dict (with birth times and birth locations)
        lin_dict = _get_lineage_dicts_one_tree(tc, tree, nodes, self.t)
        if color is None:
            # create start-color values for nodes' separate lineage tracks
            colors = [mpl.colors.to_hex(plt.cm.Set1_r(
                                        n)) for n in np.linspace(0, 0.85, 8)]
        else:
            colors = [color for _ in range(8)]
        # plot the species, either with or without phenotype-painting
        if phenotype is None:
            self._plot(lyr_num=lyr_num, land=land, size=size)
        else:
            self._plot_phenotype(phenotype, land=land, size=size)
        ax = plt.gca()
        # extract and plot the series of points for each node
        for i, node in enumerate(nodes):
            start_col = colors[i % len(colors)]
            if style == 'lineage':
                locs = np.vstack([v[1] for v in lin_dict[node].values()])
                if jitter:
                    locs = locs + np.random.normal(0, 0.01,
                                                size=locs.size).reshape(locs.shape)
                # create list of colors for plotting, using linearly interpolated
                # colors if the color argument was not provided, 
                # or else just using the solid color provided to the color argument
                color_nums = np.int8(np.linspace(0, 100, locs.shape[0]-1))
                if color is None:
                    plot_colors =[viz._calc_reshaded_color(start_col,
                                                  num) for num in color_nums]
                else:
                    plot_colors = [start_col for num in color_nums]
                # create a linear interpolation of linewidths
                linewidths = np.linspace(3, 0.85, locs.shape[0]-1)
                for n, col in enumerate(plot_colors):
                    # NOTE need to use only the first 2 values in the location
                    # data because subsequent values are used
                    # to store individuals' phenotypes and fitness
                    ax.plot(locs[n:n+2, 0], locs[n:n+2, 1], linestyle='solid',
                            marker='o', markersize=size**(1/2), color=col,
                            linewidth=linewidths[n], alpha=alpha)
            elif style == 'vector':
                # get the start and end locations
                # NOTE need to take only the first 2 values in the location
                # data because subsequent values are used to store individuals' 
                # phenotypes and fitness
                beg_loc = [*lin_dict[node].values()][-1][1][:2]
                end_loc = [*lin_dict[node].values()][0][1][:2]
                dx, dy = [end_loc[i] - beg_loc[i] for i in range(2)]
                # plot the vector
                # NOTE: SHOULD I BE FITTING A REGRESSION LINE TO THE X AND Y
                # LOCATIONS, TO GET THE 'AVERAGE' VECTOR, RATHER THAN JUST
                # PLOTTING THE VECTOR BETWEEN THE OLDEST AND CURRENT POSITIONS
                # (WHICH COULD EASILY MISREPRESENT THE OVERALL TREND BECAUSE 
                # OF CHANCE ATYPICAL LOCATIONS FOR EITHER OF THOSE TWO POSITIONS
                ax.arrow(*beg_loc, dx, dy, color=start_col,
                         width=0.05, head_width=0.4, length_includes_head=True)
            # plot the nodes' current locations and their birth locations,
            # connected by a thin black line
            node_curr_loc = node_curr_locs[i]
            node_birth_loc = [*lin_dict[node].values()][0][1]
            plt.plot([node_birth_loc[0], node_curr_loc[0]],
                     [node_birth_loc[1], node_curr_loc[1]],
                     color=start_col, linestyle=':', alpha=alpha,
                     linewidth=1.2)

        if add_roots:
            self._plot_lineage_roots(tc, tree)


    # plot the lineage for a given node and locus
    def _plot_lineage_roots(self, tc, tree, alpha=0.8, size=75):
        # get the nodes
        all_nodes = self._get_nodes()

        # get the lineage dict for all nodes
        lin_dict = _get_lineage_dicts_one_tree(tc, tree, all_nodes, self.t)

        # get the roots for all distinct lineages at this locus
        root_nodes = np.unique([[*lin_dict[n].keys()][-1] for n in all_nodes])

        # get birth locations for each root node
        root_individs = [tc.nodes[n].individual for n in root_nodes]
        root_locs = [tc.individuals[i].location for i in root_individs]

        # extract and plot the series of points for each node
        for x, y in root_locs:
            # plot the root nodes' birth locations
            plt.scatter(x, y, c='white', s=size, alpha=alpha, marker='s')
        plt.show()


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

    def _plot_pop_growth(self, expected=True, actual=True, expected_color='red',
                        actual_color='blue'):
        T = range(len(self.Nt))
        x0 = self.Nt[0] / self.K.sum()
        if expected:
            plt.plot(T, [_calc_logistic_soln(x0, self.R,
                                t) * self.K.sum() for t in T],
                     color=expected_color)
        if actual:
            plt.plot(T, self.Nt, color=actual_color)
        plt.xlabel('t')
        plt.ylabel('N(t)')

    def _plot_demographic_changes(self):
        if self._changer is None:
            print(("Species._plot_demographic_changes is not valid "
                "for species with no _SpeciesChanger object.\n"))
        else:
            self._changer._plot_dem_changes(self)


    def _plot_example_recombinant_genome(self):
        print('TODO: FIX ME!')
        return
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
        print('\t\tMAKING SPECIES %s...\n' % name, flush=True)

    #if this species should have genomes, create the genomic architecture
    if 'gen_arch' in spp_params.keys():
        # print verbose output
        if verbose:
            print('\t\t\tmaking genomic architecture...\n', flush=True)
        g_params = spp_params.gen_arch
        #make genomic_architecture
        gen_arch = _make_genomic_architecture(spp_params = spp_params,
                                                                land = land)
    else:
        gen_arch = None

    # print verbose output
    if verbose:
        print('\t\t\tmaking individuals...\n', flush=True)
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
    if spp.gen_arch is not None and not burn:
        [ind._set_z(spp.gen_arch) for ind in spp.values()]

    #make density_grid
    spp._set_dens_grids(land)

    #make movement surface, if needed
    if spp._move:
        if 'move_surf' in spp_params.movement.keys():
            if verbose:
                print(('\t\t\tmaking movement surface...\n'
                       '\t\t\t\t[can take a bit]\n'), flush=True)
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
                   '\t\t\t\t[can take a bit]\n'), flush=True)
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
                   'dispersal\n\t\t\t\t surfaces will change]\n'), flush=True)
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

