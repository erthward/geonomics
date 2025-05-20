#!/usr/bin/python
# model.py

# flake8: noqa

'''
Defines the core Model class, as well as its public and private methods
'''

#geonomics imports
from geonomics.structs.landscape import _make_landscape
from geonomics.structs.community import _make_community
from geonomics.structs.species import Species
from geonomics.structs.genome import _make_genomic_architecture
from geonomics.sim.data import (_DataCollector, _tskit_table_to_pandas,
                                _get_adhoc_sample, _format_vcf, _format_fasta)

from geonomics.sim.stats import _StatsCollector
from geonomics.utils.io import (_write_csv, _write_shapefile, _write_geojson,
                                _write_file)
from geonomics.utils._str_repr_ import _get_str_spacing
from geonomics.utils.viz import _check_display

#other imports
import numpy as np
import numpy.random as r
import random
import warnings
from copy import deepcopy
import sys, os, traceback, psutil
import matplotlib as mpl
_check_display()
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA



######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

#the Model class, which will contain the main model objects, build the
#burn-in and main functions, allow customization of the objects, manage and run
#iterations, and produce and save required stats and data
class Model:
    """
    Hierarchical data structure containing a Landscape and a Community as
    attributes, and serving as an interface for all of the core simulation
    functionality (as methods to run simulations and visualize results).
    The principal Geonomics data structure.

    Parameters
    ----------
    name : str
        A character string to be assigned as the Model's name.
    params : ParametersDict
        A ParametersDict object containing all of the parameter values needed
        instantiate the Model. (See the online documentation for a detailed
        discussion of all parameters.)
    verbose : bool, optional
        Whether or not to display information about the Model's status
        during simulations. Defaults to False.
    """

    #######################
    ### SPECIAL METHODS ###
    #######################

    def __init__(self, name, params, verbose=False):

        #get the full input params dict
        self.params = deepcopy(params)
        #get the model params (to make the following code shorter)
        m_params = self.params.model

        #set the model name (which will come from the params filename)
        if name is None:
            name = 'unnamed_model'
        self.name = name

        # get the PID of the process the model is running on
        self._pid = os.getpid()

        #set the _verbose flag
        self._verbose = verbose
        #and set the private __term_width__  and __tab_len__ attributes
        self._set_term_width()
        self._set_tab_len()

        #print verbose output
        if self._verbose:
            print('\nMAKING MODEL...\n', flush=True)

        #set the seed, if required
        self.seed = None
        if 'seed' in [*m_params]:
            self.seed = m_params.seed.num
            self._set_seeds()

        #get minimum burn-in runtime (in timesteps)
        self.burn_T = m_params.burn_T
        #and set the burn-timestep counter (burn_t) to -1 (for same reason
        #as explained for self.t)
        self.burn_t = -1
        #get main runtime (in timesteps)
        self.T = m_params.T
        #and set the timestep counter (t) to -1 (to indicate that the model
        #is unrun, and so that the first timestep will default to 0 at the 
        #beginning of the timestep
        self.t = -1

        #get the number of model iterations to run
        self.n_its = m_params.its.n_its
        #set the its list
        self.its = [*range(self.n_its)][::-1]
        #and set self.it to default to -1
        self.it = -1

        #make the land and community objects
        self.land = self._make_landscape(self._verbose)
        self.comm = self._make_community(self._verbose)

        #and make the self._never_been_run attribute False,
        #to prevent the land and community from being reset
        #on the very first iteration
        self._never_been_run = True

        #make a data._DataCollector object for the self.collecter attribute
        #if necessary
        self._data_collector = None
        if 'data' in [*m_params]:
            self._data_collector = self._make_data_collector()

        #make a stats._StatsCollector object for the self.collecter attribute
        #if necessary
        self._stats_collector = None
        if 'stats' in [*m_params]:
            self._stats_collector = self._make_stats_collector()

        #create a self.reassign_genomes attribute, which defaults to False,
        #unless any species has a genomic architecture, indicating that its
        #genomes should be reassigned after burn-in; in that case it will be 
        #reset to False as soon as the genomes are reassigned
        self.reassign_genomes = None

        # grab the rand_genarch param into a self.rand_genarch flag,
        # to determine whether or not to randomly reassign the genomes and
        # genomic architectures before each iteration
        # NOTE: only applicable if rand_comm is False; the effect will be that
        #       the burn-in will not need to be rerun,
        #       so each iteration will start
        #       with the same number and spatial distribution of individuals, 
        #       but those individuals will have different genomes each
        #       iteration (and the GenomicArchitecture object that the genomes
        #       are based upon will be different as well
        self.rand_genarch = m_params.its.rand_genarch

        #and set the orig_land and orig_comm objects, if rand_landscape and
        #rand_comm are False
        self.rand_landscape = m_params.its.rand_landscape
        self.rand_comm = m_params.its.rand_comm
        self.orig_land = None
        self.orig_comm = None
        if not self.rand_landscape:
            self.orig_land = deepcopy(self.land)
        if not self.rand_comm:
            self.orig_comm = deepcopy(self.comm)
        #and set the self.repeat_burn attribute (if this is False and
        #self.rand_comm is False, self.orig_comm will be replaced with the
        #first iteration's community as soon as it has exited burn-in and had
        #its genomes reassigned)
        self.repeat_burn = m_params.its.repeat_burn

        #set the function queues to None (they will be reset later if self.run 
        #is called, but if self.walk is called then the burn_fn_queue's None
        #status will be checked to determine whether
        #or not to run self._reset()
        self.burn_fn_queue = None
        self.main_fn_queue = None

    #define the __str__ and __repr__ special methods
    def __str__(self):
        #get a string representation of the class
        type_str = str(type(self))
        #get model name string
        name_str = "Model name:%s%s"
        name_str = name_str % (_get_str_spacing(name_str), self.name)
        #get strings for the numbers of its, burn-in timesteps, main timesteps
        n_its_str = "Number of iterations:%s%i"
        n_its_str = n_its_str % (_get_str_spacing(n_its_str), self.n_its)
        burn_t_str= "Number of burn-in timesteps (minimum):%s%i"
        burn_t_str = burn_t_str % (_get_str_spacing(burn_t_str), self.burn_T)
        main_t_str = "Number of main timesteps:%s%i"
        main_t_str = main_t_str % (_get_str_spacing(main_t_str), self.T)
        #get strings for land and comm
        comm_str = "Species:"
        spp_strs = ["%s %i: '%s'" % (_get_str_spacing(""),
            i, spp.name) for i, spp in self.comm.items()]
        spp_strs[0] = comm_str + spp_strs[0][len(comm_str):]
        comm_str = '\n'.join(spp_strs)
        land_str = "Layers:"
        lyr_strs = ["%s %i: '%s'" % (_get_str_spacing(""),
            i, lyr.name) for i, lyr in self.land.items()]
        lyr_strs[0] = land_str + lyr_strs[0][len(land_str):]
        land_str = '\n'.join(lyr_strs)
        #get string for the stats to be collected
        if self._stats_collector is not None:
            st = self._stats_collector.stats
            stat_strs=[*set([str(s) for k in st.keys() for s in st[k].keys()])]
        else:
            stat_strs = [""]
        stats_str = "Stats collected:%s{%s}"
        stats_str = stats_str % (_get_str_spacing(stats_str),
                                                    ', '.join(stat_strs))
        #get strings for the geo and gen data to be collected
        if self._data_collector is not None:
            dc = self._data_collector
            geo_data_strs = dc.geo_formats
            if dc.include_landscape:
                geo_data_strs.append(dc.rast_format)
            gen_data_strs = dc.gen_formats
        else:
            geo_data_strs = [""]
            gen_data_strs = [""]
        geo_data_str = "Geo-data collected:%s{%s}"
        geo_data_str = geo_data_str % (_get_str_spacing(geo_data_str),
                                        ', '.join(geo_data_strs))
        gen_data_str = "Gen-data collected:%s{%s}"
        gen_data_str = gen_data_str % (_get_str_spacing(gen_data_str),
                                        ', '.join(gen_data_strs))
        #join all the strings together
        tot_str = '\n'.join([type_str, name_str, land_str, comm_str,
            n_its_str, burn_t_str, main_t_str, geo_data_str, gen_data_str,
                                                            stats_str])
        return tot_str

    def __repr__(self):
        repr_str = self.__str__()
        return repr_str


    #####################
    ### OTHER METHODS ###
    #####################

        #################
        #private methods#
        #################

    #method for getting a layer's number, given a name (str) or number
    def _get_lyr_num(self, lyr_id):
        if isinstance(lyr_id, int) or lyr_id is None:
            return(lyr_id)
        elif isinstance(lyr_id, str):
            #get nums for all lyrs with matching names
            lyr_nums = [k for k, lyr  in self.land.items(
                                                    ) if lyr.name == lyr_id]
            assert len(lyr_nums) == 1, ("Expected to find a single Layer "
                "with a name matching the name provided (%s). Instead "
                "found %i.") % (lyr_id, len(lyr_nums))
            lyr_num = lyr_nums[0]
            return lyr_num
        else:
            raise ValueError(("The Layer identifier must be either a str "
                "(indicating the Layer's name) or an int (indicating "
                "its key in the Landscape dict). Instead, a %s was "
                "provided.") % str(type(lyr_id)))


    #method for getting a species' number, given a name (str) or number
    def _get_spp_num(self, spp_id):
        if isinstance(spp_id, int):
            assert spp_id in self.comm.keys(), (f"A Species with numeric "
                                                "spp_id {spp_id} does "
                                                 "not exist.")
            return(spp_id)
        elif isinstance(spp_id, str):
            #get nums for all spps with  matching names
            spp_nums = [k for k, spp in self.comm.items(
                                                    ) if spp.name == spp_id]
            assert len(spp_nums) == 1, ("Expected to find a single Species "
                "with a name matching the name provided (%s). Instead "
                "found %i.") % (spp_id, len(spp_nums))
            spp_num = spp_nums[0]
            return spp_num
        else:
            raise ValueError(("The Species identifier must be either a str "
                "(indicating the Species' name) or an int (indicating "
                "its key in the Community dict). Instead, a %s was "
                "provided.") % str(type(spp_id)))

    def _get_spp_msprime_init_status(self, spp_id):
        msprime_init_status = ('msprime' in
                        self.params['comm']['species'][spp_id]['init'])
        return msprime_init_status

    #method for getting a Trait's number, give a name (str) or number
    def _get_trt_num(self, spp, trt_id):
        if isinstance(trt_id, int) or trt_id is None:
            return(trt_id)
        elif isinstance(trt, str):
            #get nums for all traits with matching names
            trt_nums = [k for k, trt in spp.gen_arch.traits.items(
                                                    ) if trt.name == trt_id]
            assert len(trt_nums) == 1, ("Expected to find a single Trait in "
                "Species '%s' with a name matching the name provided (%s). "
                "Instead found %i.") % (spp.name, trt_id, len(trt_nums))
            trt_num = trt_nums[0]
            return trt_num
        else:
            raise ValueError(("The Trait identifier must be either a str "
                "(indicating the Trait's name) or an int (indicating "
                "its key in the Species.gen_arch.traits dict). Instead, "
                "a %s was provided.") % str(type(trt_id)))

    # private method for printing the RAM used by this model's process
    def _show_RSS(self):
        print(('\nINFO: PROCESS %i USING %0.4fGB '
               'RAM\n') % (self._pid,
                           psutil.Process(self._pid).memory_info().rss/(10**9)))
        return

    #private method for determining the width of
    #the terminal on the current system
    def _set_term_width(self):
        try:
            self.__term_width__ = int(os.popen('stty size', 'r').read().split()[1])
        except Exception:
            self.__term_width__ = 80

    #private method for determining the length of a tab on the current system
    def _set_tab_len(self):
        self.__tab_len__ = len('\t'.expandtabs())

    #method to set the iteration-counter, self.it
    def _set_it(self):
        self.it = self.its.pop()

    #method to increment self.t (i.e. the main timestep counter)
    def _set_t(self, reset=False):
        self.t += 1

    #method to reset self.t (i.e. the main timestep counter)
    def _reset_t(self):
        self.t = -1

    #method to increment self.burn_t (i.e. the burn-in timestep counter)
    def _set_burn_t(self, reset=False):
        self.burn_t += 1

    #method to reset self.burn_t (i.e. the burn-in timestep counter)
    def _reset_burn_t(self):
        self.burn_t = -1

    #method to set the self.reassign_genomes attribute 
    def _set_reassign_genomes(self):
        self.reassign_genomes = np.any(
                    [spp.gen_arch is not None for spp in self.comm.values()])

    #method to set seed (will be run when Model object is first created, if
    #called for in params)
    def _set_seeds(self):
        random.seed(self.seed)
        r.seed(self.seed)

    ###################################################
    #wrappers around functions for timesteps' fn_queues
    ###################################################

    #wrapper around Community._set_t
    def _set_comm_t(self):
        self.comm._set_t()

    #wrapper around Species._set_t
    def _set_spp_t(self, spp_idx):
        self.comm[spp_idx]._set_t()

    #wrapper around Species._set_age_stage
    def _set_age_stage(self, spp_idx):
        self.comm[spp_idx]._set_age_stage()

    #wrapper around Species._set_Nt
    def _set_Nt(self, spp_idx):
        self.comm[spp_idx]._set_Nt()

    #wrapper around Species._do_movement
    def _do_movement(self, spp_idx):
        self.comm[spp_idx]._do_movement(self.land)

    #wrapper around Species._do_pop_dynamics
    def _do_pop_dynamics(self, spp_idx):
        self.comm[spp_idx]._do_pop_dynamics(self.land)

    #wrapper around Land._make_change
    def _make_land_change(self):
        self.land._make_change(t = self.t, verbose = self._verbose)

    #wrapper around Species._make_change
    def _make_spp_change(self, spp_idx):
        self.comm[spp_idx]._make_change(verbose = self._verbose)

    #wrapper around Community._check_burned
    def _check_comm_burned(self):
        self.comm._check_burned(burn_T=self.burn_T,
                                params=self.params,
                               )

    ###################################################

    #method to wrap around landscape._make_landscape
    def _make_landscape(self, verbose=False):
        land = _make_landscape(mod=self, params=self.params,
                               verbose=verbose)
        return(land)

    #a method to reset the land as necessary
    #(recopying or regeneratiing as necessary)
    def _reset_landscape(self, rand_landscape=True):
        #deepcopy the original land if necessary 
        #(else randomly generate new land)
        if not rand_landscape:
            #verbose output
            if self._verbose:
                print(('Copying the original landscape '
                       'for iteration %i...\n\n') % self.it, flush=True)
            #deepcopy the original land
            self.land = deepcopy(self.orig_land)
            #and reset the land._changer.changes, if needed (so that they
            #point to the current land, not previous objects with possibly
            #updated attribute values)
            if self.land._changer is not None:
                #verbose output
                if self._verbose:
                    print(('Resetting the _LandscapeChanger.changes '
                                                'object...\n\n'), flush=True)
                self.land._changer._set_changes(self.land)
        else:
            #verbose output
            if self._verbose:
                print(('Creating new landscape for '
                                'iteration %i...\n\n') % self.it, flush=True)
            self.land = self._make_landscape()

    #method to wrap around Species._set_K
    def _set_K(self, spp_idx, land):
        self.comm[spp_idx]._set_K(land)

    #method to wrap around community._make_community
    def _make_community(self, verbose=False):
        comm = _make_community(self.land, self.params, burn=True,
                               verbose=verbose)
        return(comm)

    #a method to reset the community (recopying or regenerating as necessary)
    def _reset_community(self, rand_comm=True):
        if not rand_comm:
            #verbose output
            if self._verbose:
                print(('Copying the original community for '
                                'iteration %i...\n\n') % self.it, flush=True)
            self.comm = deepcopy(self.orig_comm)
            #and reset the spp._changer.changes objects for each spp,
            #if needed (so that they point to the current species,
            #not previous ones with updated attribute values)
            for spp in self.comm.values():
                if spp._changer is not None:
                    #verbose output
                    if self._verbose:
                        print(('Resetting the spp._changer.changes '
                            'object for species " %s"...\n\n') % spp.name,
                              flush=True)
                    spp._changer._set_changes(spp, self.land)
            # reset the species' GenomicArchitectures and genomes, if necessary
            if self.rand_genarch:
                for spp in self.comm.values():
                    if spp.gen_arch is not None:
                        # recreate the genomic architecture
                        if self._verbose:
                            print(('\nRecreating genomic architecture for '
                                    'species "%s"...\n\n') % spp.name,
                                  flush=True)
                        spp.gen_arch = _make_genomic_architecture(
                            spp_params=self.params['comm']['species'][spp.name],
                            land=self.land)
                        # if burn-in will not be repeated then reset the
                        # genomes and tskit tables
                        if not self.repeat_burn:
                            #verbose output
                            if self._verbose:
                                print(('\nAssigning genomes for '
                                       'species "%s"...\n\n') % spp.name,
                                      flush=True)
                            spp._set_genomes_and_tables(self.burn_T, self.T)
                            # replace gnx-simulated Individuals with
                            # msprime-simulated Individuals, if required
                            if ('msprime' in self.params['comm']['species'][
                                                            spp.name]['init']):
                                msp = self.params['comm']['species'][
                                                spp.name]['init']['msprime']
                                spp._init_msprime_pop(land=self.land,
                                                      msprime_init_params=msp,
                                                     )

        else:
            #verbose ouput
            if self._verbose:
                print(('Creating new community for '
                                    'iteration %i...\n\n') % self.it,
                      flush=True)
            self.comm = self._make_community(self._verbose)

    #method to make a data._DataCollector object for this model
    def _make_data_collector(self):
        data_collector = _DataCollector(self.name, self.params)
        return data_collector


    #method to reset the self._data_collector attribute 
    #(a data._DataCollector object)
    def _reset_data_collector(self):
        self._data_collector = self._make_data_collector()


    #method to make a stats._StatsCollector object for this model
    def _make_stats_collector(self):
        stats_collector = _StatsCollector(self.name, self.params)
        return stats_collector


    #method to reset the self._stats_collector attribute 
    #(a stats._StatsCollector object)
    def _reset_stats_collector(self):
        self._stats_collector = self._make_stats_collector()


    #method to reset all the model's objects
    #(land, community, and associated) and attributes
    def _reset(self, rand_landscape=None, rand_comm=None, rand_genarch=None,
               repeat_burn=None):
        #default to the self.rand_<_> attributes, if not otherwise provided
        if rand_landscape is None:
            rand_landscape = self.rand_landscape
        if rand_comm is None:
            rand_comm = self.rand_comm
        if rand_genarch is None:
            rand_genarch = self.rand_genarch
        if repeat_burn is None:
            repeat_burn = self.repeat_burn

        #deepcopy the original landscape and comm if necessary
        #(if told not to randomize landscape and not at the beginning
        #of the initial iteration), else randomly generate new
        if not self._never_been_run:
            self._reset_landscape(rand_landscape)
            self._reset_community(rand_comm)
        else:
            self._never_been_run = False

        #reset the self.burn_t and self.t attributes to -1
        self._reset_t()
        if repeat_burn:
            self._reset_burn_t()

        #reset the community and species t attributes
        self.comm._reset_t()
        for spp in self.comm.values():
            spp._reset_t()

        #set the self.reassign_genomes attribute
        self._set_reassign_genomes()

        #reset the self._data_collector attribute (the data._DataCollector
        #object) if necessary
        if self._data_collector is not None:
            self._reset_data_collector()

        #reset the self._stats_collector attribute (the stats._StatsCollector
        #object) if necessary
        if self._stats_collector is not None:
            self._reset_stats_collector()

        #create new main fn queue (and burn fn queue, if needed)
        if repeat_burn or self.it <= 0:
            #verbose output
            if self._verbose:
                print('Creating the burn-in function queue...\n\n', flush=True)
            self.burn_fn_queue = self._make_fn_queue(burn=True)
        #verbose output
        if self._verbose:
            print('Creating the main function queue...\n\n', flush=True)
        self.main_fn_queue = self._make_fn_queue(burn=False)


    #method to create the simulation functionality, as a function queue 
    #NOTE: (creates the list of functions that will be run by
    #self.run_burn_timestep or self.run_main_timestep,
    #depending on burn argument)
    #NOTE: because the queue is a list of functions,
    #the user could add, remove, change the order, or repeat
    #functions within the list as desired
    def _make_fn_queue(self, burn=False):
        queue = []
        #append the methods to increment the timestep
        #attributes (i.e. the timestep counters)
        if burn:
            queue.append(self._set_burn_t)
        if not burn:
            queue.append(self._set_t)
            queue.append(self._set_comm_t)
            for spp in self.comm.values():
                queue.append(lambda: self._set_spp_t(spp.idx))
        #append the set_age_stage methods to the queue
        for spp in self.comm.values():
            # NOTE: skip burn-in functions for msprime-init Species
            if not burn or not self._get_spp_msprime_init_status(spp.name):
                queue.append(lambda: self._set_age_stage(spp.idx))
        #append the do_movement_methods, if spp._move
        for spp in self.comm.values():
            # NOTE: skip burn-in functions for msprime-init Species
            if not burn or not self._get_spp_msprime_init_status(spp.name):
                if spp._move:
                    queue.append(lambda: self._do_movement(spp.idx))
        #append the do_pop_dynamics methods
        #FIXME: Consider whether the order of these needs to be specified, or
        #randomized, should people want to eventually simulate
        #multiple, interacting species
        for spp in self.comm.values():
            # NOTE: skip burn-in functions for msprime-init Species
            if not burn or not self._get_spp_msprime_init_status(spp.name):
                queue.append(lambda: self._do_pop_dynamics(spp.idx))
        #append the set_Nt methods
        # (now that this time step's population sizes are set)
        for spp in self.comm.values():
            # NOTE: skip burn-in functions for msprime-init Species
            if not burn or not self._get_spp_msprime_init_status(spp.name):
                queue.append(lambda: self._set_Nt(spp.idx))
        #add the Changer.make_change, data._DataCollector._write_data, and 
        #stats._StatsCollector._write_stats methods, if this is not the burn-in
        #and if spp and/or land have Changer objects, or if the model has 
        #_DataCollector or _StatsCollector objects (in the self._data_collector 
        #and self._stats_collector attributes)
        if not burn:
            #add land._make_change method
            if self.land._changer is not None:
                queue.append(self._make_land_change)
                #add Species._set_K methods (which will update all Species' K
                #rasters, in case landscape change has changed the Layers
                #they're based on
                for spp in self.comm.values():
                    queue.append(lambda: self._set_K(spp.idx, self.land))
            #add spp._make_change methods
            for spp in self.comm.values():
                if spp._changer is not None:
                    queue.append(lambda: self._make_spp_change(spp.idx))
            #add self._write_data method
            if self._data_collector is not None:
                queue.append(self._write_data)
            #add self.calc_stats method
            if self._stats_collector is not None:
                queue.append(self.calc_stats)

        #add the burn-in function if need be
        if burn:
            queue.append(self._check_comm_burned)
        return(queue)


    #method to generate timestep_verbose_output
    def _print_timestep_info(self, mode):
        verbose_msg = '%s:\tit=%i:\tt=%i\n' % (mode, self.it,
                                [self.burn_t if mode == 'burn' else self.t][0])
        spps_submsgs = ''
        for spp in self.comm.values():
            if len(spp.Nt) > 0:
                Nt_val = spp.Nt[:].pop()
            else:
                Nt_val = np.nan
            if len(spp.n_births) > 0:
                n_births_val = spp.n_births[:].pop()
            else:
                n_births_val = np.nan
            if len(spp.n_deaths) > 0:
                n_deaths_val = spp.n_deaths[:].pop()
            else:
                n_deaths_val = np.nan
            spacer = ' ' * (30 - len(spp.name))
            spp_submsg = (f"\tspecies: {spp.name}{spacer}N={Nt_val}"
                          f"\t(births={n_births_val}\tdeaths={n_deaths_val})\n")
            spps_submsgs += spp_submsg
        verbose_msg = verbose_msg + spps_submsgs
        print(verbose_msg)
        print('\t' + '.' * (self.__term_width__ - self.__tab_len__), flush=True)


    #method to run the burn-in or main function
    #queue (as indicated by mode argument)
    def _do_timestep(self, mode):
        #do a burn-in step
        if mode == 'burn':
            for fn in self.burn_fn_queue:
                if True not in [spp.extinct for spp in self.comm.values()]:
                    fn()
                else:
                    break
            if self._verbose:
                self._print_timestep_info(mode)
            # if the burn-in is complete, reassign the genomes if needed,
            # set the tskit.TableCollection tables if needed,
            # and then set self.comm.burned = True
            if np.all([spp.burned for spp in self.comm.values()]):
                if self.reassign_genomes:
                    for spp in self.comm.values():
                        if spp.gen_arch is not None:
                            #verbose output
                            if self._verbose:
                                print(('\nAssigning genomes for '
                                       'species "%s"...\n\n') % spp.name,
                                      flush=True)
                            spp._set_genomes_and_tables(self.burn_T, self.T)
                            # replace gnx-simulated Individuals with
                            # msprime-simulated Individuals, if required
                            if self._get_spp_msprime_init_status(spp.name):
                                msp = self.params['comm']['species'][
                                                spp.name]['init']['msprime']
                                spp._init_msprime_pop(land=self.land,
                                                      msprime_init_params=msp,
                                                     )
                    #and then set the reassign_genomes attribute to False, so
                    #that they won't get reassigned again during this iteration
                    self.reassign_genomes = False
                #mark the community as burned in
                self.comm.burned = True
                #verbose output
                if self._verbose:
                    print('Burn-in complete.\n\n', flush=True)
        #or do a main step
        elif mode == 'main':
            for fn in self.main_fn_queue:
                if  True not in [spp.extinct for spp in self.comm.values()]:
                    fn()
                else:
                    break
            # set each species' flag that indicates whether or not the
            # TableCollection is sorted and simplified to False
            for spp in self.comm.values():
                if (spp.gen_arch is not None
                    and spp.gen_arch.use_tskit
                    and spp._tc_sorted_and_simplified):
                    spp._tc_sorted_and_simplified = False

            if self._verbose:
                self._print_timestep_info(mode)
            # sort and simplify tskit tables, if needed
            for spp in self.comm.values():
                    if (spp.gen_arch is not None
                        and spp.gen_arch.use_tskit
                        and (self.t + 1) % spp.gen_arch.tskit_simp_interval == 0
                        and self.t !=-1):
                        if self._verbose:
                            print(('\n\nnow sorting and simplifying '
                                   "tskit tables for Species: '%s'") % spp.name)
                            print("\tNUMBER EDGES BEFORE SIMPLIFICATION:",
                                  self.comm[0]._tc.edges.num_rows)
                            print("\tNUMBER INDIVIDS BEFORE SIMPLIFICATION:",
                                  self.comm[0]._tc.individuals.num_rows)
                        spp._sort_and_simplify_table_collection()
                        if self._verbose:
                            print("\n\tNUMBER EDGES AFTER SIMPLIFICATION: ",
                                  self.comm[0]._tc.edges.num_rows)
                            print("\tNUMBER INDIVIDS AFTER SIMPLIFICATION: ",
                                  self.comm[0]._tc.individuals.num_rows, '\n',
                                  flush=True)

        #then check if any species are extinct and
        #return the correpsonding boolean
        extinct = np.any([spp.extinct for spp in self.comm.values()])
        #verbose output
        if extinct and self._verbose:
            print(('XXXX     Species %s went extinct. '
                'Iteration %i aborting.\n\n') % (' & '.join(
                ['"' + spp.name + '"' for spp in self.comm.values(
                                                ) if spp.extinct]), self.it),
                 flush=True)

        return(extinct)


    #method to set the model's next iteration
    def _set_next_iteration(self, core=None):
        #update the iteration numbers
        self._set_it()
        #verbose output
        if self._verbose:
            print('~' * self.__term_width__ + '\n\n')
            print('Setting up iteration %i...\n\n' % self.it, flush=True)

        #reset the model (including the burn-in,
        #if self.repeat_burn or if this is the first iteration) 
        self._reset(rand_landscape = self.rand_landscape,
                    rand_comm = self.rand_comm,
                    rand_genarch=self.rand_genarch,
                    repeat_burn = self.repeat_burn)


    #method for running the model's next iteration
    def _do_next_iteration(self):

        #set the next iteration to be run
        self._set_next_iteration()

        #loop over the burn-in timesteps running the burn-in
        #queue (if copied spp isn't already burned in)
        if (self.rand_comm
            or (not self.rand_comm and self.repeat_burn)
            or self.it == 0):
            #verbose output
            if self._verbose:
                print('Running burn-in, iteration %i...\n\n' % self.it,
                      flush=True)
            #until all species have spp.burned == True
            while not np.all([spp.burned for spp in self.comm.values()]):
                #run a burn-in timestep
                extinct = self._do_timestep(mode = 'burn')
                #and end the iteration early if any species is extinct
                if extinct:
                    break

            #overwrite self.orig_comm with the burned-in community,
            #if necessary
            if not self.rand_comm and not self.repeat_burn:
                #verbose output
                if self._verbose:
                    print('Saving burned-in community...\n\n', flush=True)
                self.orig_comm = deepcopy(self.comm)

        #verbose output
        if self._verbose:
            print('Running main model, iteration %i...\n\n' % self.it,
                  flush=True)

        #check if any species is starting out extinct (which would happen if it
        #went extinct in the burn-in), and if so, then skip main mode
        if np.any([spp.extinct for spp in self.comm.values()]):
            if self._verbose:
                print(("WARNING: At least one Species went extinct during "
                    "the burn-in. Cannot run main phase for "
                    "iteration %i.\n\n") % self.it, flush=True)
            return

        #loop over the timesteps, running the run_main function repeatedly
        for t in range(self.T):
            #run a main timestep
            extinct = self._do_timestep('main')
            #and end the iteration early if any species is extinct 
            if extinct:
                break


        ##################
        # public methods #
        ##################

    # method to run the overall model
    def run(self, verbose=False):
        """
        Run a Model.

        Run the Model for the number of iterations that was stipulated in the
        parameters file used to create it (which is saved as the `n_its`
        attribute). During each iteration, the model will burn in for at least
        the stipulated number of burn-in timesteps (the `burn_T` attribute),
        then will run for the stipulated number of main timesteps (the 'T'
        attribute). If stats and/or data are going to be collected (also
        stipulated by the parameters file used to create the Model), then the
        output files for each iteration will be saved in a separate
        subdirectory, containing further subdirectories for each Species.

        Parameters
        ----------
        verbose : bool, optional
            Whether or not to run the Model should provide written output.
            If True, formatted messages will be printed to STDOUT at each
            timestep (which could be piped to a log file, if desired).

        Returns
        -------
        out : None
            Returns no output. Writes information to STDOUT if verbose is
            True. Writes data and statistics to file if so parameterized.

        Examples
        --------
        We can create and run the default model (as long as there is no
        "GNX_params_<...>.py" file in the current working directory
        before `gnx.make_parameters_file` is called).

        >>> gnx.make_parameters_file()
        >>> mod = gnx.make_model()
        NOTE: Using the following file, in the current working directory to create the Model object:
                GNX_params_21-01-2020_17:22:08.py
        >>> mod.run(verbose = True)
        ###############################################################################################
        Running model "GNX_params_21-01-2020_17:22:08"...
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Setting up iteration 0...
        Creating the burn-in function queue...
        Creating the main function queue...
        Running burn-in, iteration 0...
        burn:   it=0:   t=0
                species: spp_0                         N=250    (births=35      deaths=214)
        .......................................................................................
        burn:   it=0:   t=1
                species: spp_0                         N=250    (births=35      deaths=214)
        .......................................................................................
        burn:   it=0:   t=2
                species: spp_0                         N=250    (births=35      deaths=214)
        .......................................................................................
        .
        .
        .
        """
        # TODO: Add some handler for farming instances out to
        # different nodes, if available?

        # set the self._verbose flag
        self._verbose = verbose

        # verbose output
        if self._verbose:
            print('\n\n' + '#' * self.__term_width__ + '\n\n')
            print('Running model "%s"...\n\n' % self.name, flush=True)

        # loop over all the iterations
        while len(self.its) > 0:
            # do the next iteration
            try:
                self._do_next_iteration()
            except Exception as e:
                msg = ('XXXX\tAn error occurred during iteration %i, '
                       'timestep %i.\n \tPLEASE COPY AND REPORT THE FOLLOWING,'
                       ' to help us debug Geonomics:\n')
                if self.comm.burned:
                    msg = msg % (self.it, self.t)
                else:
                    msg = msg % (self.it, self.burn_t)
                print(msg)
                print('<>' * int(self.__term_width__/2))
                print('Error message:\n\t%s\n\n' % e)
                traceback.print_exc(file=sys.stdout)
                print()
                print('<>' * int(self.__term_width__/2) + '\n\n', flush=True)

        # verbose output
        if self._verbose:
            print('\n\nModel "%s" is complete.\n' % self.name)
            print('#' * self.__term_width__, flush=True)

        # set the _verbose flag back to False
        self._verbose = False


    # method to run the model interactively from the command line; named 'walk'
    # to succinctly differentiate it from 'run'
    def walk(self, T=1, mode='main', verbose=True, animate=False):
        """
        Walk through a Model (i.e. run it for a certain number of timesteps).

        This function will run a Model for a specific number of timesteps,
        either in 'burn' or 'main' mode. It is designed to help users to try
        out, explore, and introspect their Models in an interactive Python
        session. (For best results, we strongly recommend this be done with
        the `iPython <https://ipython.org>` shell.)

        Parameters
        ----------
        T : int, optional
            The number of timesteps for which to run the Model. (If in 'main'
            mode, the Model will run for precisely this many timesteps. If in
            'burn' mode, it will run either until stationarity statistics
            indicate that it has burned in, or for this many timesteps,
            whichever occurs first.) Defaults to a single timestep.
            single step.
        mode : {'burn', 'main'}, optional
            The mode in which to run the Model.
            If 'main', it will run in main mode (i.e. all functionalities will
            be used. (Note that for a Model to run in main mode it must first
            have been burned in.)
            If 'burn', it will run in burn-in mode. This entails the following:
            1. Genomes will not yet be assigned, so genomic phenomena such
            as crossing-over and natural selection will not occur;
            2. Any scheduled landscape and demographic changes will not occur;
            3. Stationarity statistics will be checked at the end of each
            timestep, to determine if the Model has burned in (in which
            case the model will stop running regardless of timestep,
            genomes will be randomly assigned to all individuals, and
            the Model will be deemed ready to run in main mode).
        verbose : bool, optional
            Whether or not to run the Model should provide written output.
            If True, formatted messages will be printed to STDOUT at each
            timestep. Defaults to True for `Model.walk`.
        animate : {tuple of ints, bool}, optional
            If a length-2 tuple of integers is provided, the Species indicated
            by the first number will be plotted on the Landscape Layer
            indicated by the second number at each timestep, in a dynamically
            updating plot.
            If a length-3 tuple of integers in provided, the Species indicated
            by the first number will be plotted on the Landscape Layer
            indicated by the second number at each timestep, colored by the
            phenotypes of the trait indicated by the third number.
            If just True, the first Species (index 0) will be plotted
            on top of a transparent stack of all Landscape Layers.
            (Note that this will slow down the execution of a model a bit,
            because the plot will pause for 0.1 seconds after each timestep.)

        Returns
        -------
        out : None
            Returns no output. Writes information to STDOUT if verbose is
            True. Writes data and statistics to file if so parameterized.

        Raises
        ------
        ValueError
            If the user attempts to run the Model in 'main' mode before it has
            been burned in.

        Examples
        --------
        We can create and walk the default model (as long as there is no
        "GNX_params_<...>.py" file in the current working directory
        before `gnx.make_parameters_file` is called).

        >>> gnx.make_parameters_file()
        >>> mod = gnx.make_model()
        NOTE: Using the following file, in the current working directory to create the Model object:
                GNX_params_21-01-2020_17:22:08.py
        >>> #run the burn-in until it is complete
        >>> mod.walk(T=1000, mode='burn')
        ###############################################################################################
        Running model "GNX_params_21-01-2020_17:22:08"...
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Setting up iteration 0...
        Creating the burn-in function queue...
        Creating the main function queue...
        Running burn-in, iteration 0...
        burn:   it=0:   t=0
                species: spp_0                         N=250    (births=35      deaths=214)
        .......................................................................................
        burn:   it=0:   t=1
                species: spp_0                         N=250    (births=35      deaths=214)
        .......................................................................................
        burn:   it=0:   t=2
                species: spp_0                         N=250    (births=35      deaths=214)
        .......................................................................................
        .
        .
        .
        Assigning genomes for species "spp_0"...
        Burn-in complete.
        >>> #now run in main mode for 50 timesteps
        >>> mod.walk(T=50, mode='main')
         main:   it=1:  t=0
                 species: spp_0                         N=131    (births=35      deaths=214)
         .......................................................................................
         main:   it=1:  t=1
                 species: spp_0                         N=129    (births=35      deaths=214)
         .......................................................................................
         main:   it=1:  t=2
                 species: spp_0                         N=133    (births=35      deaths=214)
         .......................................................................................
         .
         .
         .
             """
        # validate T, and convert to int if a float fed (such as 1e6,
        # to run burn-in to completion)
        assert isinstance(T, int) or isinstance(T, float), ("'T' must be "
                                                            "a numeric data "
                                                            "type.")
        if isinstance(T, float):
            T = int(T)

        # throw an error and quit if mode == 'main' but the community is
        # not burned in
        if mode == 'main' and not self.comm.burned:
            raise ValueError(("The Model.walk method cannot be run in 'main' "
                              "mode if the Model's Community has not yet "
                              "been burned in (i.e. if Model.comm.burned "
                              "is False)."))

        # set the model's _verbose flag
        old_verbose = deepcopy(self._verbose)
        self._verbose = verbose

        # if verbose, add a space below the command line prompt,
        # for readability
        if self._verbose:
            print('\n')

        # coerce animate to tuple if list was provided
        if isinstance(animate, list):
            animate = tuple(animate)
        # start animated plot, if animate == True
        if animate not in (False, None):
            if (isinstance(animate, tuple)
               and sum([isinstance(i, int) for i in animate]) == len(animate)):
                if len(animate) == 2:
                    points = self.plot(animate[0], animate[1], animate=True)
                elif len(animate) == 3:
                    points = self.plot_phenotype(animate[0], animate[1],
                                                 animate[2], animate=True)
            elif animate is True:
                points = self.plot(spp=0, animate=True)
            plt.ion()
            plt.draw()
            plt.pause(0.05)

        # run the model for the stipulated number of timesteps
        for t in range(T):
            # exit if burn-in is complete
            if mode == 'burn' and self.comm.burned:
                # verbose output
                break
            # reset mode, if mode is 'burn' and there is no mod.burn_fn_queue
            if mode == 'burn' and self.burn_fn_queue is None:
                # verbose output
                if self._verbose:
                    print(('No mod.burn_fn_queue was found. '
                          'Running mod.reset()...\n\n'), flush=True)
                self._reset()
            extinct = self._do_timestep(mode=mode)
            # continue the animated plot, if animate == True
            if animate not in (False, None):
                if (isinstance(animate, tuple)
                   and sum([isinstance(i,
                                       int) for i in animate]) == len(
                                                                    animate)):
                    points.remove()
                    if len(animate) == 2:
                        points = self.plot(animate[0], animate[1],
                                           animate=True)
                    elif len(animate) == 3:
                        points = self.plot_phenotype(animate[0], animate[1],
                                                     animate[2], animate=True)
                    plt.ion()
                    plt.draw()
                    plt.pause(0.1)
                elif animate in [True]:
                    points.remove()
                    points = self.plot(spp=0, animate=True)
                    plt.ion()
                    plt.draw()
                    plt.pause(0.1)

            # end the iteration early if any species is extinct
            if extinct:
                break
        # reset self._verbose to False
        self._verbose = old_verbose


    # helper method to burn a model in if it hasn't been burned in already
    # to succinctly differentiate it from 'run'
    def burn(self):
        """
        Helper function to burn a model in, if it hasn't been burned in already.

        Returns
        -------
        out : None
            Returns no output. Alters Model in place by running burn-in phase.

        """
        if self.comm.burned:
            print("\nModel has already been burned in.\n")
        else:
            while not self.comm.burned:
                self.walk(10, 'burn')
            print("\nModel has now been burned in.\n")


    # method to use the self._data_collector object to sample and write data
    def _write_data(self):
        self._data_collector._write_data(self.comm, self.land, self.it)


    #method to use the self._stats_collector object to sample and write data
    def calc_stats(self):
        self._stats_collector._calc_stats(self.comm, self.t, self.it)


    ##########
    #plotting#
    ##########

    #wrapper around Species._plot and Landscape._plot
    #TODO: allow spp to be a list of spp ids, to plot each spp in a
    #different color with a legend!
    def plot(self, spp=None, lyr=None, hide_land=False, individs=None,
             text=False, color='black', edge_color='face', text_color='black',
             cbar=True, size=25, text_size=9, land_cmap=None, pt_cmap=None,
             alpha=False, zoom_width=None, x=None, y=None, vmin=None,
             vmax=None, ticks=False, mask_rast=None, animate=False):
        """
        Create a basic plot of the Model in its current state.

        Plot the current state of the model, with any number of overlain
        Landscape Layers, and with or without points representing the
        Individuals in a Species.

        Parameters
        ----------
        spp : {int, str}, optional, default: None
            A reference to the Species whose Individuals should be scattered
            on the plot. Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string). If None, will cause only the Landscape to be plotted.

        lyr : {int, str}, optional, default: None
            A reference to the Layer whose raster should be plotted. Can be
            either the Layer's index number (i.e. its integer key in the
            Landscape dict), or its name (as a character string). Defaults
            to None, which will cause all Layers to be plotted as an overlay
            of transparent rasters, each with a different colormap.

        hide_land : bool, optional, default: False
            If True, the Landscape will be 'hidden', i.e. the plot will only
            show a scatter of the plotted Individuals. Defaults to False.

        individs : iterable collection of ints, optional, default: None
            If provided, indicates the integer indices of the only Individuals
            to be plotted. If None, all Individuals will be plotted.

        text : bool, optional, default: False
            If True, each Individual's index number will be displayed next to
            it. Can be useful for model introspection.

        color : valid mpl.plt color value, optional, default: 'black'
            Face color for the points in the Individual scatter. Passed to the
            *c* argument of matplotlib.pyplot.scatter.

        edge_color : valid mpl.plt color value, optional, default: 'face'
            Edge color for the points in the Individual scatter. If 'face',
            will always match the face color (i.e. the color provided to
            *color*). Passed to the *edgecolor* argument of
            matplotlib.pyplot.scatter.

        text_color : valid mpl.plt color value, default 'black'
            Color for the plotted text. (Will only be used if the *text*
            argument is True.) Passed to the *color* argument of
            matplotlib.pyplot.text.

        cbar : bool, optional, default: True
            If True, a colorbar will be included, depicting the mapping of
            color onto the Landscape's environmental values.

        size : scalar or array_like, optional, default: 25
            Size of the scatter points. Passed to the *s* argument
            of matplotlib.pyplot.scatter.

        text_size : {size in points, valid string}, optional, default: 9
            Text size. (Will only be used if *text* is True.) Can be expressed
            as a numeric size in points, or as any string that is a valid size
            argument for matplotlib.text.Text.

        land_cmap : {valid string, None}, optional, default: None
            Colormap to use for plotting Layer rasters. If None, will default
            to the colormap automatically assigned to the Layer by Geonomics
            based on the Layer's index number in the Landscape dict. Can be
            passed any string that references a colormap in
            matplotlib.pyplot.cm.

        pt_cmap : {valid string, None}, optional, default: None
            Colormap to use for plotting points. (Will only be used if the
            Individual scatter is being plotted, and if the value provided to
            *color* is an ordered collection of values to which to map the
            Individuals' colors.) If None, no colormap will be used.
            Can be passed any string that references a colormap in
            matplotlib.pyplot.cm, as well as an
            matplotlib.colors.LinearSegmentedColormap instance.

        alpha : {scalar, bool, None}, optional, default: None
            The transparency level of the points. If scalar is passed, it is
            passed on to the *alpha* value fed to matplotlib.pyplot.scatter.
            If True, alpha value will be set to 0.6. If False or None, points
            will be opaque.

        zoom_width : {scalar, None}, optional, default: None
            Width, in raster cells, of the window to which to zoom to resulting
            plot. If None, plot shows full Landscape.

        x,y : {scalar, None}, default: None
            The x and y coordinates of the center of the resulting plot.
            Only used if *zoom_width* is not None.

        vmin, vmax : {scalar}, optional, default: None
            Used to normalize the luminance data that determine how colors
            are colormapped onto Individuals' point colors. Passed to the
            *vmin* and *vmax* arguments of matplotlib.pyplot.scatter.

        ticks : bool, optional, default: False
            If True, x- and y-axis ticks will be added to the plot.

        mask_rast : array_like, optional default: None
            An array to use to mask values in the Layer(s) being plotted.
            All np.nan values in the array will be plotted as light gray,
            instead of the color that would otherwise map to their
            environmental values in their raster(s).
            The array must be of the same x,y dimensions as the Layer(s)
            being plotted.

        animate : bool, optional, default: False
            If True, this function will return a list of the Individuals'
            points, as a matplotlib.collections.PathCollection. This is used
            internally by the Model.walk method, to create animated simulation
            plots.

        Returns
        -------
        {None, :class:`matplotlib.collections.PathCollection`}
            Returns no output (unless *animate* is passed True, in which case
            a PathCollection of the Individuals' points is returned, for
            internal use by Model.walk, to create animated simulations)

        Notes
        -----
        - For more detail on plotting parameters, see the documentation for
          matplotlib, matplotlib.pyplot.pcolormesh, matplotlib.pyplot.scatter,
          and matplotlib.pyplot.text.

        """
        #get the lyr num
        lyr_num = self._get_lyr_num(lyr)
        #if no spp provided, then call Landscape._plot
        if (spp is None or
            (spp is not None and len(self.comm[self._get_spp_num(spp)])) == 0):
            self.land._plot(lyr_num=lyr_num, cbar=cbar, cmap=land_cmap, x=x,
                            y=y, zoom_width=zoom_width, vmin=vmin, vmax=vmax,
                            ticks=ticks, mask_rast=mask_rast)
            points = None
        #or else plot the spp
        else:
            #get the spp
            spp = self.comm[self._get_spp_num(spp)]
            #feed args into spp._plot
            points = spp._plot(lyr_num=lyr_num, land=self.land,
                               hide_land=hide_land, individs=individs,
                               text=text, color=color, edge_color=edge_color,
                               text_color=text_color, cbar=cbar, size=size,
                               text_size=text_size,
                               land_cmap=land_cmap, pt_cmap=pt_cmap,
                               alpha=alpha, zoom_width=zoom_width, x=x, y=y,
                               vmin=vmin, vmax=vmax, ticks=ticks,
                               mask_rast=mask_rast, animate=animate)
            #add spp name
            #plt.suptitle(spp.name)
        return points


    #wrapper around Species._plot_density
    def plot_density(self, spp=0, normalize=False, individs=None,
            text=False, color='black', edge_color='face',
            text_color='black', size=25, text_size = 9,
            alpha=0.5, zoom_width=None, x=None, y=None, ticks=None,
            mask_rast=None):
        """
        Plot the population-density raster for a given Species.

        Plot the population density of a Species, as estimated by the Species'
        _DensityGridStack.

        Parameters
        ----------
        spp : {int, str}, optional, default: 0
            A reference to the Species whose Individuals should be scattered
            on the plot. Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string).

        normalize : bool, optional, default: False
            If True, the population-density values of the output raster will be
            normalized to 0 <= density <= 1. Otherwise, actual estimates
            (individuals per cell) are returned.

        individs : iterable collection of ints, optional, default: None
            If provided, indicates the integer indices of the only Individuals
            to be plotted. If None, all Individuals will be plotted.

        text : bool, optional, default: False
            If True, each Individual's index number will be displayed next to
            it. Can be useful for model introspection.

        color : valid mpl.plt color value, optional, default: 'black'
            Face color for the points in the Individual scatter. Passed to the
            *c* argument of matplotlib.pyplot.scatter.

        edge_color : valid mpl.plt color value, optional, default: 'face'
            Edge color for the points in the Individual scatter. If 'face',
            will always match the face color (i.e. the color provided to
            *color*). Passed to the *edgecolor* argument of
            matplotlib.pyplot.scatter.

        text_color : valid mpl.plt color value, default 'black'
            Color for the plotted text. (Will only be used if the *text*
            argument is True.) Passed to the *color* argument of
            matplotlib.pyplot.text.

        size : scalar or array_like, optional, default: 25
            Size of the scatter points. Passed to the *s* argument
            of matplotlib.pyplot.scatter.

        text_size : {size in points, valid string}, optional, default: 9
            Text size. (Will only be used if *text* is True.) Can be expressed
            as a numeric size in points, or as any string that is a valid size
            argument for matplotlib.text.Text.

        alpha : {scalar, bool, None}, optional, default: None
            The transparency level of the points. If scalar is passed, it is
            passed on to the *alpha* value fed to matplotlib.pyplot.scatter.
            If True, alpha value will be set to 0.6. If False or None, points
            will be opaque.

        zoom_width : {scalar, None}, optional, default: None
            Width, in raster cells, of the window to which to zoom to resulting
            plot. If None, plot shows full Landscape.

        x,y : {scalar, None}, default: None
            The x and y coordinates of the center of the resulting plot.
            Only used if *zoom_width* is not None.

        ticks : bool, optional, default: False
            If True, x- and y-axis ticks will be added to the plot.

        mask_rast : array_like, optional default: None
            An array to use to mask values in the Layer(s) being plotted.
            All np.nan values in the array will be plotted as light gray,
            instead of the color that would otherwise map to their
            environmental values in their raster(s).
            The array must be of the same x,y dimensions as the Layer(s)
            being plotted.

        Returns
        -------
        None
            Returns no output

        Notes
        -----
        - For more detail on plotting parameters, see the documentation for
          matplotlib, matplotlib.pyplot.pcolormesh, matplotlib.pyplot.scatter,
          and matplotlib.pyplot.text.

        """

        #get the spp
        spp = self.comm[self._get_spp_num(spp)]
        #feed args into spp._plot_density
        spp._plot_density(land = self.land, normalize=normalize,
            individs=individs, text=text, color=color, edge_color=edge_color,
            text_color=text_color, size=size, text_size=text_size, alpha=alpha,
            zoom_width=zoom_width, x=x, y=y, ticks=ticks, mask_rast=mask_rast)
        #add spp name
        #plt.suptitle(spp.name)


    #wrapper around Species._plot_genotype
    def plot_genotype(self, spp=0, locus=0, lyr=None, by_dominance=False,
                      individs=None, text=False, edge_color='black',
                      text_color='black', cbar=True, size=25, text_size = 9,
                      alpha=1, zoom_width=None, x=None, y=None, ticks=None,
                      mask_rast=None):
        """
        Plot a Species' Individuals, colored by their genotypes for a locus

        Create a scatter plot of the Individuals in a Species,
        on top of any Landscape Layer, and color the points by
        their genotypes for the chosen locus.

        Parameters
        ----------
        spp : {int, str}, optional, default: 0
            A reference to the Species whose Individuals should be scattered
            on the plot. Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string).

        locus : int, optional, default: 0
            The index of the locus for which Individuals' genotypes should be
            depicted. Must be 0 <= val <= (genome length - 1)

        lyr : {int, str}, optional, default: None
            A reference to the Layer whose raster should be plotted. Can be
            either the Layer's index number (i.e. its integer key in the
            Landscape dict), or its name (as a character string). Defaults
            to None, which will cause all Layers to be plotted as an overlay
            of transparent rasters, each with a different colormap.

        by_dominance : bool, optional, default: False
            If False, the Individuals will be colored by their actual genotypes
            (i.e. a 0|0 homozygote will be depicted as 0, a 1|1 homozygote as
            1, and a 0|1 heterozygote as 0.5). If True, then the dominance
            value for the chosen locus will be used to the determine the
            'effective genotype', which will instead be displayed as an
            Individual's color. In this case, 0|1 heterozygotes
            will be be depicted as 0.5's for loci with codominance (i.e.
            dominance values of 0), but
            will be depicted as 1's for loci with 1-allele dominance (i.e.
            dominance values of 1). The dominance values for all loci in the
            genome are set in a Species' GenomicArchitecture object. For more
            detail, see the "Genomes, GenomicArchitecture, and Traits"
            subsection of the "Data Structures" section of the online
            documentation.

        individs : iterable collection of ints, optional, default: None
            If provided, indicates the integer indices of the only Individuals
            to be plotted. If None, all Individuals will be plotted.

        text : bool, optional, default: False
            If True, each Individual's index number will be displayed next to
            it. Can be useful for model introspection.

        edge_color : valid mpl.plt color value, optional, default: 'face'
            Edge color for the points in the Individual scatter. If 'face',
            will always match the face color (i.e. the color provided to
            *color*). Passed to the *edgecolor* argument of
            matplotlib.pyplot.scatter.

        text_color : valid mpl.plt color value, default 'black'
            Color for the plotted text. (Will only be used if the *text*
            argument is True.) Passed to the *color* argument of
            matplotlib.pyplot.text.

        cbar : bool, optional, default: True
            If True, a colorbar will be included, depicting the mapping of
            color onto the Landscape's environmental values.

        size : scalar or array_like, optional, default: 25
            Size of the scatter points. Passed to the *s* argument
            of matplotlib.pyplot.scatter.

        text_size : {size in points, valid string}, optional, default: 9
            Text size. (Will only be used if *text* is True.) Can be expressed
            as a numeric size in points, or as any string that is a valid size
            argument for matplotlib.text.Text.

        alpha : {scalar, bool, None}, optional, default: None
            The transparency level of the points. If scalar is passed, it is
            passed on to the *alpha* value fed to matplotlib.pyplot.scatter.
            If True, alpha value will be set to 0.6. If False or None, points
            will be opaque.

        zoom_width : {scalar, None}, optional, default: None
            Width, in raster cells, of the window to which to zoom to resulting
            plot. If None, plot shows full Landscape.

        x,y : {scalar, None}, default: None
            The x and y coordinates of the center of the resulting plot.
            Only used if *zoom_width* is not None.

        ticks : bool, optional, default: False
            If True, x- and y-axis ticks will be added to the plot.

        mask_rast : array_like, optional default: None
            An array to use to mask values in the Layer(s) being plotted.
            All np.nan values in the array will be plotted as light gray,
            instead of the color that would otherwise map to their
            environmental values in their raster(s).
            The array must be of the same x,y dimensions as the Layer(s)
            being plotted.

        Returns
        -------
        None
            Returns no output

        Notes
        -----
        - For more detail on plotting parameters, see the documentation for
          matplotlib, matplotlib.pyplot.pcolormesh, matplotlib.pyplot.scatter,
          and matplotlib.pyplot.text.

        """
        #get the lyr num
        lyr_num = self._get_lyr_num(lyr)
        #get the spp
        spp = self.comm[self._get_spp_num(spp)]
        #feed args into spp._plot_genotype
        spp._plot_genotype(locus=locus, lyr_num=lyr_num, land=self.land,
                           individs=individs, text=text, size=size,
                           text_size=text_size,
                           edge_color=edge_color, text_color=text_color,
                           cbar=cbar, alpha=alpha,
                           by_dominance=by_dominance,
                           zoom_width=zoom_width, x=x, y=y, ticks=ticks,
                           mask_rast=None)
        #add spp name
        #plt.suptitle(spp.name)


    # wrapper around Species._plot_phenotype
    # for a given trait
    def plot_phenotype(self, spp=0, trt=0, lyr=None, individs=None,
                       text=False,  edge_color='black', text_color='black',
                       cbar=True, size=25, text_size=9, alpha=1,
                       zoom_width=None, x=None, y=None, ticks=None,
                       mask_rast=None, animate=False):
        """
        Plot a Species' Individuals, colored by their phenotypes for a Trait

        Create a scatter plot of the Individuals in a Species,
        on top of any Landscape Layer, and color the points by
        their phenotypes for the chosen Trait.

        Parameters
        ----------
        spp : {int, str}, optional, default: 0
            A reference to the Species whose Individuals should be scattered
            on the plot. Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string).

        trt : {int, str}, optional, default: 0
            A reference to the Trait whose phenotypes should be used to color
            Individuals in the plot. Can be either the Trait's index number
            (i.e. its integer key in the GenomicArchitecture's traits dict),
            or its name (as a character string).

        lyr : {int, str}, optional, default: None
            A reference to the Layer whose raster should be plotted. Can be
            either the Layer's index number (i.e. its integer key in the
            Landscape dict), or its name (as a character string). Defaults
            to None, which will cause only the Layer associated with the chosen
            Trait to be plotted.

        individs : iterable collection of ints, optional, default: None
            If provided, indicates the integer indices of the only Individuals
            to be plotted. If None, all Individuals will be plotted.

        text : bool, optional, default: False
            If True, each Individual's index number will be displayed next to
            it. Can be useful for model introspection.

        edge_color : valid mpl.plt color value, optional, default: 'face'
            Edge color for the points in the Individual scatter. If 'face',
            will always match the face color (i.e. the color provided to
            *color*). Passed to the *edgecolor* argument of
            matplotlib.pyplot.scatter.

        text_color : valid mpl.plt color value, default 'black'
            Color for the plotted text. (Will only be used if the *text*
            argument is True.) Passed to the *color* argument of
            matplotlib.pyplot.text.

        cbar : bool, optional, default: True
            If True, a colorbar will be included, depicting the mapping of
            color onto the Landscape's environmental values.

        size : scalar or array_like, optional, default: 25
            Size of the scatter points. Passed to the *s* argument
            of matplotlib.pyplot.scatter.

        text_size : {size in points, valid string}, optional, default: 9
            Text size. (Will only be used if *text* is True.) Can be expressed
            as a numeric size in points, or as any string that is a valid size
            argument for matplotlib.text.Text.

        alpha : {scalar, bool, None}, optional, default: None
            The transparency level of the points. If scalar is passed, it is
            passed on to the *alpha* value fed to matplotlib.pyplot.scatter.
            If True, alpha value will be set to 0.6. If False or None, points
            will be opaque.

        zoom_width : {scalar, None}, optional, default: None
            Width, in raster cells, of the window to which to zoom to resulting
            plot. If None, plot shows full Landscape.

        x,y : {scalar, None}, default: None
            The x and y coordinates of the center of the resulting plot.
            Only used if *zoom_width* is not None.

        ticks : bool, optional, default: False
            If True, x- and y-axis ticks will be added to the plot.

        mask_rast : array_like, optional default: None
            An array to use to mask values in the Layer(s) being plotted.
            All np.nan values in the array will be plotted as light gray,
            instead of the color that would otherwise map to their
            environmental values in their raster(s).
            The array must be of the same x,y dimensions as the Layer(s)
            being plotted.

        animate : bool, optional, default: False
            If True, this function will return a list of the Individuals'
            points, as a matplotlib.collections.PathCollection. This is used
            internally by the Model.walk method, to create animated simulation
            plots.

        Returns
        -------
        {None, :class:`matplotlib.collections.PathCollection`}
            Returns no output (unless *animate* is passed True, in which case
            a PathCollection of the Individuals' points is returned, for
            internal use by Model.walk, to create animated simulations)

        Notes
        -----
        - For more detail on plotting parameters, see the documentation for
          matplotlib, matplotlib.pyplot.pcolormesh, matplotlib.pyplot.scatter,
          and matplotlib.pyplot.text.

        """

        # get the lyr num
        lyr_num = self._get_lyr_num(lyr)
        # get the spp
        spp = self.comm[self._get_spp_num(spp)]
        # return messages if species does not have genomes or traits
        if spp.gen_arch is None:
            print(("Model.plot_phenotype is not valid for Species "
                "without genomes.\n"), flush=True)
            return
        elif spp.gen_arch.traits is None:
            print(("Model.plot_phenotype is not valid for Species "
                "without traits.\n"), flush=True)
            return
        # get the trt_num
        trt_num = self._get_trt_num(spp, trt)
        # trt_num can't be None for plot_phenotype
        assert trt_num is not None, ("None is not a valid value for the "
            "'trait' arguemnt.")
        # feed args into spp._plot_phenotype
        points = spp._plot_phenotype(trait=trt, lyr_num=lyr_num,
                                     land=self.land, individs=individs,
                                     text=text, size=size, text_size=text_size,
                                     edge_color=edge_color,
                                     text_color=text_color, cbar=cbar,
                                     alpha=alpha, zoom_width=zoom_width,
                                     x=x, y=y, ticks=ticks,
                                     mask_rast=mask_rast, animate=animate)
        #add spp name
        #plt.suptitle(spp.name)
        return points


    #wrapper around Species._plot_fitness
    def plot_fitness(self, spp=0, trt=None, lyr=None, individs=None,
                     text=False, phenotype_text=False, fitness_text=False,
                     edge_color='black', text_color='black',
                     phenotype_text_color='black', fitness_text_color='black',
                     cbar=True, fitness_cbar=True, size=100, text_size=9,
                     fit_cmap='gray', alpha=1, zoom_width=None, x=None, y=None,
                     ticks=None, mask_rast=None):
        """
        Plot a Species' Individuals, colored by fitness

        Create a scatter plot of the Individuals in a Species,
        on top of any Landscape Layer, and color the points by
        their fitnesses. Fitness values can be either overall fitness
        (i.e. calculated across all traits), or can be fitnesses for
        a certain Trait.

        Parameters
        ----------
        spp : {int, str}, optional, default: 0
            A reference to the Species whose Individuals should be scattered
            on the plot. Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string). If None, will cause only the Landscape
            to be plotted.

        trt : {int, str}, optional, default: None
            A reference to the Trait for which fitness should be calculated
            and used to color Individuals in the plot. Can be either the
            Trait's index number (i.e. its integer key in the
            GenomicArchitecture's traits dict), or its name (as a character
            string). If None, Individuals' overall fitnesses will be
            calculated instead.

        lyr : {int, str}, optional, default: None
            A reference to the Layer whose raster should be plotted. Can be
            either the Layer's index number (i.e. its integer key in the
            Landscape dict), or its name (as a character string). Defaults
            to None, which will either cause only the Layer associated with the
            chosen Trait to be plotted, if *trt* is not None; or else
            will cause all Layers to be plotted as an overlay of transparent
            rasters, each with a different colormap, if *trt* is None.

        individs : iterable collection of ints, optional, default: None
            If provided, indicates the integer indices of the only Individuals
            to be plotted. If None, all Individuals will be plotted.

        text : bool, optional, default: False
            If True, each Individual's index number will be displayed next to
            it. Can be useful for model introspection.

        phenotype_text : bool, optional, default: False
            If True, each Individual's phenotype will be displayed next to
            it. Can be useful for model introspection.

        fitness_text : bool, optional, default: False
            If True, each Individual's calculated fitness value will be
            displayed next to it. Can be useful for model introspection.

        edge_color : valid mpl.plt color value, optional, default: 'face'
            Edge color for the points in the Individual scatter. If 'face',
            will always match the face color (i.e. the color provided to
            *color*). Passed to the *edgecolor* argument of
            matplotlib.pyplot.scatter.

        text_color : valid mpl.plt color value, default 'black'
            Color for the plotted text. (Will only be used if the *text*
            argument is True.) Passed to the *color* argument of
            matplotlib.pyplot.text.

        phenotype_text_color, fitness_text_color : valid mpl.plt color value,
        default 'black'
            Color for the plotted phenotype-value or fitness-value text.
            (Will only be used if the *phenotype_text* or *fitness_text*
            arguments are True.) Passed to the *color* argument
            of matplotlib.pyplot.text.

        cbar : bool, optional, default: True
            If True, a colorbar will be included, depicting the mapping of
            color onto the Landscape's environmental values.

        fitness_cbar : bool, optional, default: True
            If True, a colorbar will be included, depicting the mapping of
            color onto the Landscape's environmental values.

        size : scalar or array_like, optional, default: 25
            Size of the scatter points. Passed to the *s* argument
            of matplotlib.pyplot.scatter.

        text_size : {size in points, valid string}, optional, default: 9
            Text size. (Will only be used if *text* is True.) Can be expressed
            as a numeric size in points, or as any string that is a valid size
            argument for matplotlib.text.Text.

        fit_cmap : {valid string, None}, optional, default: 'gray'
            Colormap to use for color Individuals by their fitness values.
            If no value passed, will defaul to grayscale. Can be passed any
            string that references a colormap in matplotlib.pyplot.cm.

        alpha : {scalar, bool, None}, optional, default: None
            The transparency level of the points. If scalar is passed, it is
            passed on to the *alpha* value fed to matplotlib.pyplot.scatter.
            If True, alpha value will be set to 0.6. If False or None, points
            will be opaque.

        zoom_width : {scalar, None}, optional, default: None
            Width, in raster cells, of the window to which to zoom to resulting
            plot. If None, plot shows full Landscape.

        x,y : {scalar, None}, default: None
            The x and y coordinates of the center of the resulting plot.
            Only used if *zoom_width* is not None.

        ticks : bool, optional, default: False
            If True, x- and y-axis ticks will be added to the plot.

        mask_rast : array_like, optional default: None
            An array to use to mask values in the Layer(s) being plotted.
            All np.nan values in the array will be plotted as light gray,
            instead of the color that would otherwise map to their
            environmental values in their raster(s).
            The array must be of the same x,y dimensions as the Layer(s)
            being plotted.

        Returns
        -------
        None
            Returns no output

        Notes
        -----
        - For more detail on plotting parameters, see the documentation for
          matplotlib, matplotlib.pyplot.pcolormesh, matplotlib.pyplot.scatter,
          and matplotlib.pyplot.text.

        """

        #get the lyr num
        lyr_num = self._get_lyr_num(lyr)
        #get the spp
        spp = self.comm[self._get_spp_num(spp)]
        #return messages if species does not have genomes or traits
        if spp.gen_arch is None:
            print(("Model.plot_fitness is not valid for Species "
                "without genomes.\n"), flush=True)
            return
        elif spp.gen_arch.traits is None:
            print(("Model.plot_fitness is not valid for Species "
                "without traits.\n"), flush=True)
            return
        #get the trt_num, which CAN be None for plot_fitness
        trt_num = self._get_trt_num(spp, trt)
        #feed args into spp._plot_fitness
        spp._plot_fitness(trt_num=trt_num, lyr_num=lyr_num, land = self.land,
                          individs=individs, text=text,
                          phenotype_text=phenotype_text,
                          phenotype_text_color=phenotype_text_color,
                          fitness_text=fitness_text,
                          fitness_text_color=fitness_text_color, size=size,
                          text_size=text_size, edge_color=edge_color,
                          text_color=text_color, fit_cmap=fit_cmap, cbar=cbar,
                          fitness_cbar=fitness_cbar, alpha=alpha,
                          zoom_width=zoom_width, x=x, y=y, ticks=ticks,
                          mask_rast=mask_rast)
        #add spp name
        #plt.suptitle(spp.name)


    def plot_genetic_PCA(self, n_pcs=3, spp=0, lyr=0,
                        individs=None, text=False,
                        edge_color='black', text_color='black', cbar=True,
                        size=25, text_size=9,
                        alpha=False, zoom_width=None, x=None, y=None,
                        ticks=False, hide_land=False, mask_rast=None):
        """
        Create a scatter plot of the Individuals in a Species,
        on top of any Landscape Layer, and color the points by
        the first 3 PCS in a genetic PCA (Principal Components Analysis).

        Parameters
        ----------
        spp : {int, str}, optional, default: 0
            A reference to the Species whose Individuals should be scattered
            on the plot. Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string). If None, will cause only the Landscape
            to be plotted.

        n_pcs : int, optional, default: 3
            The number of top PCs to use for the plot. If 1 is used, all three
            red, green, and blue (RGB) channels will be scaled the same,
            creating a grayscale mapping of the top PC. If 2 are used, the red
            and blue channels only will be used, creating a visualization of
            the top 2 PCs that is colorblind-friendly. If 3 are used, the top 3
            PCs will be fed into all RGB channels, creating a visualization
            with the most information but which is not colorblind-friendly.

        lyr : {int, str}, optional, default: None
            A reference to the Layer whose raster should be plotted. Can be
            either the Layer's index number (i.e. its integer key in the
            Landscape dict), or its name (as a character string). Defaults
            to None, which will either cause only the Layer associated with the
            chosen Trait to be plotted, if *trt* is not None; or else
            will cause all Layers to be plotted as an overlay of transparent
            rasters, each with a different colormap, if *trt* is None.

        individs : iterable collection of ints, optional, default: None
            If provided, indicates the integer indices of the only Individuals
            to be plotted. If None, all Individuals will be plotted.

        text : bool, optional, default: False
            If True, each Individual's index number will be displayed next to
            it. Can be useful for model introspection.

        edge_color : valid mpl.plt color value, optional, default: 'face'
            Edge color for the points in the Individual scatter. If 'face',
            will always match the face color (i.e. the color provided to
            *color*). Passed to the *edgecolor* argument of
            matplotlib.pyplot.scatter.

        text_color : valid mpl.plt color value, default 'black'
            Color for the plotted text. (Will only be used if the *text*
            argument is True.) Passed to the *color* argument of
            matplotlib.pyplot.text.

        cbar : bool, optional, default: True
            If True, a colorbar will be included, depicting the mapping of
            color onto the Landscape's environmental values.

        size : scalar or array_like, optional, default: 25
            Size of the scatter points. Passed to the *s* argument
            of matplotlib.pyplot.scatter.

        text_size : {size in points, valid string}, optional, default: 9
            Text size. (Will only be used if *text* is True.) Can be expressed
            as a numeric size in points, or as any string that is a valid size
            argument for matplotlib.text.Text.

        alpha : {scalar, bool, None}, optional, default: None
            The transparency level of the points. If scalar is passed, it is
            passed on to the *alpha* value fed to matplotlib.pyplot.scatter.
            If True, alpha value will be set to 0.6. If False or None, points
            will be opaque.

        zoom_width : {scalar, None}, optional, default: None
            Width, in raster cells, of the window to which to zoom to resulting
            plot. If None, plot shows full Landscape.

        x,y : {scalar, None}, default: None
            The x and y coordinates of the center of the resulting plot.
            Only used if *zoom_width* is not None.

        ticks : bool, optional, default: False
            If True, x- and y-axis ticks will be added to the plot.

        hide_land : bool, optional, default: False
            If True, Landscape Layer raster will not be plotted beneath the
            scatterplot.

        mask_rast : array_like, optional default: None
            An array to use to mask values in the Layer(s) being plotted.
            All np.nan values in the array will be plotted as light gray,
            instead of the color that would otherwise map to their
            environmental values in their raster(s).
            The array must be of the same x,y dimensions as the Layer(s)
            being plotted.

        Returns
        -------
        None
            Returns no output

        Notes
        -----
        - For more detail on plotting parameters, see the documentation for
          matplotlib, matplotlib.pyplot.pcolormesh, matplotlib.pyplot.scatter,
          and matplotlib.pyplot.text.

        """
        #get the lyr num
        lyr_num = self._get_lyr_num(lyr)
        #get the spp
        spp = self.comm[self._get_spp_num(spp)]

        # get array of resulting genomic data (i.e. 'speciome'),
        # genotypes meaned by individual
        if spp.gen_arch.use_tskit:
            speciome = spp._get_genotypes(biallelic=False)
        else:
            speciome = np.mean(np.stack([i.g for i in spp.values()]), axis=2)
        assert (0 < n_pcs <= 3 and
                isinstance(n_pcs, int)), "'n_pcs' must be an integer between 0 and 3"
        # run PCA on speciome
        pca = PCA(n_components=n_pcs)
        PCs = pca.fit_transform(speciome)
        # normalize the PC results
        norm_PCs = (PCs - np.min(PCs,
                                 axis=0)) / (np.max(PCs,
                                                    axis=0) - np.min(PCs,
                                                                     axis=0))
        # use top PCs to get normalized values for R, G, & B colors
        # grayscale if n_pcs==1
        # red-blue if n_pcs==2
        # RGB if n_pcs==3
        if n_pcs == 1:
            PC_colors = np.hstack([norm_PCs * 255 for _ in range(3)])
        if n_pcs == 2:
            PC_colors = np.vstack([norm_PCs[:, 0] * 255,
                                   np.zeros(norm_PCs.shape[0]),
                                   norm_PCs[:, 1] * 255,
                                  ]).T
        if n_pcs == 3:
            PC_colors = norm_PCs * 255

        # scatter all individuals on top of landscape, colored by top PCs
        # (number determined by n_pcs)
        xs = self.comm[0]._get_x()
        ys = self.comm[0]._get_y()

        # subset xs, ys, and PC_colors, if individs requested
        if individs is not None:
            indices = [i in individs for i in spp]
            PC_colors = PC_colors[indices, :]
            xs = xs[indices]
            ys = ys[indices]

        # plot it
        self.plot(spp=spp.idx, lyr=lyr_num, color=PC_colors/255.0,
                  individs=individs, text=text,
                  edge_color=edge_color, text_color=text_color, cbar=cbar,
                  size=size, text_size=text_size, alpha=alpha,
                  zoom_width=zoom_width, x=x, y=y, ticks=ticks,
                  hide_land=hide_land, mask_rast=mask_rast)


    #wrapper around Species._plot_allele_frequencies
    def plot_allele_frequencies(self, spp=0, color='red'):
        """
        Plot a Species' current allele frequencies for all alleles

        For the Species indicated, create a plot showing, across all loci,
        the current allele frequencies (in the indicated color) versus the starting
        frequencies (in dashed black).

        Parameters
        ----------
        spp : {int, str}, optional, default: 0
            A reference to the Species whose allele frequencies will be
            plotted. Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string).

        color: str, optional, default: 'red'
            The color to use for plotting the current allele frequencies.

        Returns
        -------
        None
            Returns no output.
        """
        #get the spp
        spp = self.comm[self._get_spp_num(spp)]
        #call the fn
        spp._plot_allele_frequencies(color=color)


    #wrapper around Species._plot_hist_fitness
    def plot_hist_fitness(self, spp=0):
        """
        Plot a Species

        For the Species indicated, create a histogram of Individuals' current
        fitness values.

        Parameters
        ----------
        spp : {int, str}, optional, default: 0
            A reference to the Species whose Individuals' fitness will be
            plotted. Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string).

        Returns
        -------
        None
            Returns no output.
        """
        #get the spp
        spp = self.comm[self._get_spp_num(spp)]
        #call the fn
        spp._plot_hist_fitness()


    def _plot_direction_surface(self, surf_type, spp, style, x=None, y=None,
                                zoom_width=None, scale_fact=4.5, color='black',
                                cbar = True, ticks=None, cmap='plasma',
                                mask_rast=None, head_width=0.24,
                                head_length=0.32):
        '''
        wrapper around Species._plot_direction_surface
        '''
        #get the spp
        spp = self.comm[self._get_spp_num(spp)]
        #call the fn
        spp._plot_direction_surface(land=self.land, surf_type=surf_type,
            style=style, x=x, y=y, zoom_width=zoom_width,
            scale_fact=scale_fact, color=color, cbar=cbar, ticks=ticks,
            cmap=cmap, mask_rast=mask_rast, head_width=head_width,
                                    head_length=head_length)


    #wrapper around Species._plot_direction_surface for _move_surf
    def plot_movement_surface(self, spp=0, style='hist',  color='black', cbar=True,
                              ticks=None, cmap='plasma', zoom_width=None,
                              x=None, y=None, scale_fact=4.5, mask_rast=None,
                              head_width=0.24, head_length=0.32):
        """
        Create any of 4 types of plots to summarize a Species' MovementSurface

        Plot a Species' MovemementSurface with either circularized histograms,
        a circle of directional draws, or a vector on each cell, to visualize
        the von Mises distributions of each cell on the MovementSurface.
        Or plot a plain histogram of the von Mises distribution at a given
        cell.

        Parameters
        ----------
        spp : {int, str}, optional, default: 0
            A reference to the Species whose MovementSurface should be plotted.
            Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string).

        style : str, optional, default: 'hist'
            Determines the style of plot produced. There are four options:
                'hist':
                    Plot a classic histogram approximating the von Mises
                    distribution at the cell indicated by position x,y.
                    (Note that *x* and *y* must be provided for this style.)
                'chist':
                    Plot a circular histogram approximating the von Mises
                    distribution at the cell indicated by position x,y;
                    plot will be drawn inside the chosen cell on the
                    MovementSurface raster.
                'cdraws':
                    Plot points on the unit circle, whose locations were
                    drawn from the von Mises distribution at the cell indicated
                    by position x,y; plot will be drawn inside the chosen cell
                    on the MovementSurface raster.
                'vect':
                    Inside each cell of the MovementSurface raster,
                    plot the mean direction vector of directions drawn
                    from that cell's von Mises distribution.

        color : valid mpl.plt color value, optional, default: 'black'
            Color to use for the plot elements. Passed to the
            *c* argument of matplotlib.pyplot.hist, matplotlib.pyplot.plot,
            matplotlib.pyplot.scatter, or matplotlib.pyplot.arrow,
            depending on *style*.

        cbar : bool, optional, default: True
            If True, a colorbar will be included, depicting the mapping of
            color onto the MovementSurface's environmental values.

        cmap : {valid string, None}, optional, default: 'plasma'
            Colormap to use for plotting the MovementSurface. If None,
            will default to 'plasma', a common colormap used to plot resistance
            surfaces.Can be passed any string that references a colormap in
            matplotlib.pyplot.cm.

        zoom_width : {scalar, None}, optional, default: None
            Width, in raster cells, of the window to which to zoom to resulting
            plot. If None, plot shows full Landscape.

        x,y : {scalar, None}, default: None
            The x and y coordinates of the center of the resulting plot.
            Only used if *zoom_width* is not None.

        mask_rast : array_like, optional default: None
            An array to use to mask values in the Layer(s) being plotted.
            All np.nan values in the array will be plotted as light gray,
            instead of the color that would otherwise map to their
            environmental values in their raster(s).
            The array must be of the same x,y dimensions as the Layer(s)
            being plotted.

        Returns
        -------
        None
            Returns no output

        Notes
        -----
        - Not valid for Species without a DisersalSurface
        - For more detail on plotting parameters, see the documentation for
          matplotlib, matplotlib.pyplot.pcolormesh, matplotlib.pyplot.scatter,
          and matplotlib.pyplot.text.

        """

        self._plot_direction_surface(surf_type='move', spp=spp, style=style,
            x=x, y=y, zoom_width=zoom_width, scale_fact=scale_fact,
            color=color, cbar=cbar, ticks=ticks, cmap=cmap,
            mask_rast=mask_rast, head_width=head_width, head_length=head_length)


    #wrapper around Species._plot_direciton_surface for _disp_surf
    def plot_dispersal_surface(self, spp=0, style='hist', x=None, y=None,
                               zoom_width=None, scale_fact=4.5, color='black',
                               cbar = True, ticks=None, cmap='plasma',
                               mask_rast=None, head_width=0.24,
                               head_length=0.32):
        """
        Create any of 4 types of plots to summarize a Species' DispersalSurface

        Plot a Species' DispersalSurface with either circularized histograms,
        a circle of directional draws, or a vector on each cell, to visualize
        the von Mises distributions of each cell on the MovementSurface.
        Or plot a plain histogram of the von Mises distribution at a given
        cell.

        Parameters
        ----------
        spp : {int, str}, optional, default: 0
            A reference to the Species whose DispersalSurface should be
            plotted. Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string).

        style : str, optional, default: 'hist'
            Determines the style of plot produced. There are four options:
                'hist':
                    Plot a classic histogram approximating the von Mises
                    distribution at the cell indicated by position x,y.
                    (Note that *x* and *y* must be provided for this style.)
                'chist':
                    Plot a circular histogram approximating the von Mises
                    distribution at the cell indicated by position x,y;
                    plot will be drawn inside the chosen cell on the
                    DispersalSurface raster.
                'cdraws':
                    Plot points on the unit circle, whose locations were
                    drawn from the von Mises distribution at the cell indicated
                    by position x,y; plot will be drawn inside the chosen cell
                    on the DispersalSurface raster.
                'vect':
                    Inside each cell of the DispersalSurface raster,
                    plot the mean direction vector of directions drawn
                    from that cell's von Mises distribution.

        color : valid mpl.plt color value, optional, default: 'black'
            Color to use for the plot elements. Passed to the
            *c* argument of matplotlib.pyplot.hist, matplotlib.pyplot.plot,
            matplotlib.pyplot.scatter, or matplotlib.pyplot.arrow,
            depending on *style*.

        cbar : bool, optional, default: True
            If True, a colorbar will be included, depicting the mapping of
            color onto the MovementSurface's environmental values.

        cmap : {valid string, None}, optional, default: 'plasma'
            Colormap to use for plotting the DispersalSurface. If None,
            will default to 'plasma', a common colormap used to plot resistance
            surfaces.Can be passed any string that references a colormap in
            matplotlib.pyplot.cm.

        zoom_width : {scalar, None}, optional, default: None
            Width, in raster cells, of the window to which to zoom to resulting
            plot. If None, plot shows full Landscape.

        x,y : {scalar, None}, default: None
            The x and y coordinates of the center of the resulting plot.
            Only used if *zoom_width* is not None.

        mask_rast : array_like, optional default: None
            An array to use to mask values in the Layer(s) being plotted.
            All np.nan values in the array will be plotted as light gray,
            instead of the color that would otherwise map to their
            environmental values in their raster(s).
            The array must be of the same x,y dimensions as the Layer(s)
            being plotted.

        Returns
        -------
        None
            Returns no output

        Notes
        -----
        - Not valid for Species without a DisersalSurface
        - For more detail on plotting parameters, see the documentation for
          matplotlib, matplotlib.pyplot.pcolormesh, matplotlib.pyplot.scatter,
          and matplotlib.pyplot.text.

        """
        self._plot_direction_surface(surf_type='move', spp=spp, style=style,
            x=x, y=y, zoom_width=zoom_width, scale_fact=scale_fact,
            color=color, cbar=cbar, ticks=ticks, cmap=cmap,
            mask_rast=mask_rast, head_width=head_width, head_length=head_length)


    def _plot_movement_or_dispersal(self,
                                    what,
                                    spp=0,
                                    lyr=0,
                                    n_timesteps=1,
                                    n_individs=None,
                                    distance_distr_param1=None,
                                    distance_distr_param2=None,
                                    distance_distr='lognormal',
                                    color='black',
                                    alpha=None,
                                    land_cmap='plasma',
                                    size=10,
                                    ticks=False,
                                    color_by_individ=True,
                                    increasing_linewidth=True,
                                    include_start_points=False,
                                   ):
        '''
        wrapper around Species._plot_movement_or_dispersal
        '''
        assert what in ['movement', 'dispersal']
        #get the spp and the lyr-num
        spp = self.comm[self._get_spp_num(spp)]
        lyr_num = self._get_lyr_num(lyr)
        #call the fn
        spp._plot_movement_or_dispersal(what=what,
                                        land=self.land,
                                        lyr=lyr_num,
                                        n_timesteps=n_timesteps,
                                        n_individs=n_individs,
                                        distance_distr_param1=distance_distr_param1,
                                        distance_distr_param2=distance_distr_param2,
                                        distance_distr=distance_distr,
                                        color=color,
                                        alpha=alpha,
                                        land_cmap=land_cmap,
                                        size=size,
                                        ticks=ticks,
                                        color_by_individ=color_by_individ,
                                        increasing_linewidth=increasing_linewidth,
                                        include_start_points=include_start_points,
                                       )


    def plot_movement(self,
                      spp=0,
                      lyr=0,
                      n_timesteps=1,
                      n_individs=None,
                      distance_distr_param1=None,
                      distance_distr_param2=None,
                      distance_distr='lognormal',
                      color='black',
                      alpha=None,
                      land_cmap='plasma',
                      size=10,
                      ticks=False,
                      color_by_individ=True,
                      increasing_linewidth=True,
                      include_start_points=False,
                     ):
        '''
        plot example movement tracks for a number of
        individuals of the indicated species and using the indicated parameter
        values. note that directions will be drawn using whatever directional
        behavior is parameterized for the given species (i.e., from a surface,
        or simply from a universally parameterized von mises distribution).

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which movement should be plotted.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string). if None, the first species will be plotted to be plotted.

        lyr : {int, str}, optional, default: None
            a reference to the layer whose raster should be plotted. can be
            either the layer's index number (i.e. its integer key in the
            landscape dict), or its name (as a character string). defaults
            to None, in which case the first layer will be plotted.

        n_timesteps : int
            a number of timesteps (i.e. individual movement steps) to be plotted
            in each individual's movement track.

        n_individs : {int, None}, default: None
            the number of individuals whose movement should be
            plotted. defaults to None, in which case all individuals will be
            plotted.

        distance_distr_param1,
        distance_distr_param2 : float, optional default: None
            these are the parameters defining the distributions from which the
            direction and distance values for individual movements are drawn. if
            passed, their values will be used to determine the movement behavior
            that is plotted (even if a parametersdict was passed to *params*,
            because these values will override any values therein). for more
            details, see the 'parameters' section of the online documentation.

        distance_distr: str, optional default: 'lognormal'
            this indicates the distribution to use for drawing movement
            distancs. defaults to 'lognormal', but can also be set to 'wald' or
            'levy'.

        color : valid mpl.plt color value, optional, default: 'black'
            face color for the lines in the movement-track plot.

        alpha : valid mpl.plt alpha value, optional, default: None
            transparency value for the lines in the movement-track plot.

        size : scalar or array_like, optional, default: 10
            size of the scatter points.

        land_cmap : {valid string, None}, optional, default: 'plasma'
            colormap to use for plotting layer rasters. If None, will default
            to the colormap automatically assigned to the layer by geonomics
            based on the layer's index number in the landscape dict. can be
            passed any string that references a colormap in matplotlib.pyplot.cm.

        ticks : bool, optional, default: False
            If True, x- and y-axis ticks will be added to the plot.

        color_by_individ : bool, optional, default: True
            If True, each plotted individual's movement track will have its own
            distinct color, making it easier to follow the track and distinguish it
            from others in the plot.

        increasing_linewidth : bool, optional, default: True
            If True, each individual's movement track will get wider with each
            successive movement step, making it easier to distinguish
            directionality in the movement track.

        include_start_points : bool, optional, default: False
            If True, a scatter of points will be plotted to signify the starting
            points of individuals' movement tracks.
        '''
        self._plot_movement_or_dispersal(what='movement',
                                         spp=spp,
                                         lyr=lyr,
                                         n_timesteps=n_timesteps,
                                         n_individs=n_individs,
                                         distance_distr_param1=distance_distr_param1,
                                         distance_distr_param2=distance_distr_param2,
                                         distance_distr=distance_distr,
                                         color=color,
                                         alpha=alpha,
                                         land_cmap=land_cmap,
                                         size=size,
                                         ticks=ticks,
                                         color_by_individ=color_by_individ,
                                         increasing_linewidth=increasing_linewidth,
                                         include_start_points=include_start_points,
                                        )


    def plot_dispersal(self,
                       spp=0,
                       lyr=0,
                       n_dispersals=1,
                       n_individs=None,
                       distance_distr_param1=None,
                       distance_distr_param2=None,
                       distance_distr='lognormal',
                       color='black',
                       alpha=None,
                       land_cmap='plasma',
                       size=10,
                       ticks=False,
                       color_by_individ=True,
                       increasing_linewidth=True,
                       include_start_points=False,
                      ):
        '''
        plot example dispersal events for a number of
        individuals of the indicated species and using the indicated parameter
        values. note that directions will be drawn using whatever directional
        behavior is parameterized for the given species (i.e., from a surface,
        or simply from a universally parameterized von mises distribution).

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which dispersal should be plotted.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string). If None, the first species will be plotted to be plotted.

        lyr : {int, str}, optional, default: None
            a reference to the layer whose raster should be plotted. can be
            either the layer's index number (i.e. its integer key in the
            landscape dict), or its name (as a character string). defaults
            to None, in which case the first layer will be plotted.

        n_dispersals : int
            the number of independent dispersals to plot for each individual.

        n_individs : {int, None}, default: None
            the number of individuals whose dispersal should be
            plotted. defaults to None, in which case all individuals will be
            plotted.

        distance_distr_param1,
        distance_distr_param2 : float, optional default: None
            these are the parameters defining the distributions from which the
            direction and distance values for individual dispersals are drawn. if
            passed, their values will be used to determine the dispersal behavior
            that is plotted (even if a parametersdict was passed to *params*,
            because these values will override any values therein). for more
            details, see the 'parameters' section of the online documentation.

        distance_distr: str, optional default: 'lognormal'
            this indicates the distribution to use for drawing dispersal
            distancs. defaults to 'lognormal', but can also be set to 'wald' or
            'levy'.

       color : valid mpl.plt color value, optional, default: 'black'
            face color for the dipseral lines.

        alpha : valid mpl.plt alpha value, optional, default: None
            transparency value for the dispersal lines.

        size : scalar or array_like, optional, default: 10
            size of the scatter points.

        land_cmap : {valid string, None}, optional, default: 'plasma'
            colormap to use for plotting layer rasters. If None, will default
            to the colormap automatically assigned to the layer by geonomics
            based on the layer's index number in the landscape dict. can be
            passed any string that references a colormap in matplotlib.pyplot.cm.

        ticks : bool, optional, default: False
            If True, x- and y-axis ticks will be added to the plot.

        color_by_individ : bool, optional, default: True
            If True, each dispersal line for each plotted individual will have its own
            distinct color, making it easier to follow the track and distinguish it
            from others in the plot.

        include_start_points : bool, optional, default: False
            If True, a scatter of points will be plotted to signify the starting
            points from which each individual is dispersing.
        '''
        self._plot_movement_or_dispersal(what='dispersal',
                                         spp=spp,
                                         lyr=lyr,
                                         n_timesteps=n_dispersals,
                                         n_individs=n_individs,
                                         distance_distr_param1=distance_distr_param1,
                                         distance_distr_param2=distance_distr_param2,
                                         distance_distr=distance_distr,
                                         color=color,
                                         alpha=alpha,
                                         land_cmap=land_cmap,
                                         size=size,
                                         ticks=ticks,
                                         color_by_individ=color_by_individ,
                                         increasing_linewidth=False,
                                         include_start_points=include_start_points,
                                        )


    #wrapper around species._plot_demographic_pyramid
    def plot_demographic_pyramid(self, spp=0):
        """
        plot a demographic pyramid for the chosen species

        plot a paired, horizontal bar plot of the chosen species'
        count of individuals of each age.

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species whose demographic pyramid will be
            plotted. can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        returns
        -------
        None
            returns no output
        """
        #get the spp
        spp = self.comm[self._get_spp_num(spp)]
        #call the fn
        spp._plot_demographic_pyramid()


    #wrapper around species._plot_pop_growth
    def plot_pop_growth(self, spp=0, expected=True, actual=True,
                        expected_color='red', actual_color='blue'):
        """
        plot the chosen species' population over time

        plot a line plot of the species' population over model time (in blue),
        as well as the expected population-size trend (based on summing the
        species' carrying-capacity raster and using that total expected
        population size to solve the the logistic growth equation; in red).

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species whose population growth will be plotted.
            can be either the species' index number (i.e. its integer key
            in the community dict), or its name (as a character string).
        expected: bool, default: True
            whether or not to plot the line showing the expected (under
            the logistic growth model) trajectory of populaion growth, in red.
        actual: bool, default: True
            whether or not to plot the line showing the actual
            trajectory of population growth, in blue.

        returns
        -------
        None
            returns no output
        """
        #get the spp
        spp = self.comm[self._get_spp_num(spp)]
        #call the fn
        spp._plot_pop_growth(expected=expected, actual=actual,
                             expected_color=expected_color,
                             actual_color=actual_color)


    #wrapper around species._plot_example_recombinant_genome
    def plot_example_recombinant_genome(self, spp=0):
        """
        plot an example recombinant genome for the chosen species

        plot a vertical image depicting an example recombinant genome for the
        chosen species. each gamete (i.e. vertical half) is colored using a
        different colormap, and for each one alternating bands of color
        indicate breakpoints where recombination occurred between the two
        parental chromatids (or haploid chromatid sets, for multi-chromosome
        models). the resulting plot can be useful for introspecting or
        exploring a model, because each plot will be unique, and will be a
        result of the array of interlocus recombination rates defined in the
        species' genomicarchitecture.

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which an example recombinant genome
            should be plotted. can be either the species' index number
            (i.e. its integer key in the community dict), or its name
            (as a character string).

        returns
        -------
        None
            returns no output

        notes
        -----
        - not valid for species without genomes.
        """
        spp = self.comm[self._get_spp_num(spp)]
        spp._plot_example_recombinant_genome()



    ##########
    #analysis#
    ##########

    #wrapper around species._run_cca and other gea methods
    def run_gea(self, method='cca', spp=0, trt=0,
                plot=True, plot_sd=True, scale=3, sd=3):
        """
        runs one of a number of available gea models, returning the results and
        optionally plotting them.

        parameters
        ----------
        method : str, optional, default: "cca"
            a string indicating which type of gea model to run.
            note: currently, only cca (default) is available.

        spp : {int, str}, optional, default: 0
            a reference to the species for which the gea should be run.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        trt : {int, str}, optional, default: 0
            a reference to the trait for which the gea should be run.
            can be either the trait's index number (i.e. its integer key in the
            genomicarchitecture's traits dict), or its name (as a character
            string).

        plot : bool, optional, default: True
            whether or not the gea model should be plotted.

        plot_sd: bool, optional, default: True
            whether or not a standard-deviation ellipse should be added to the
            gea's plot. (only used if plot==True and method=="cca".)

        scale : {int, float}, optional, default: 3
            the scaling factor to use when plotting the gea results. (only used
            if plot==True.)

        sd : {int, float}, optional, default: 3
            the number of standard deviations to use when plotting the
            standard-deviation ellipse. (only used if plot==True and
            method=="cca".)

        returns
        -------
        dict
            returns a dict containing the labeled (as keys) output data
            structures (as values) of the chosen method of gea.

        """
        # list of valid methods
        methods = ["cca"]

        #get the spp
        spp = self.comm[self._get_spp_num(spp)]

        # get the trait number
        trt_num = self._get_trt_num(spp, trt)

        #feed args into the appropriate gea function
        if method.lower() == 'cca':
            results = spp._run_cca(trt_num=trt_num, scale=scale,
                                   plot=plot, plot_sd=plot_sd, sd=sd)
        else:
            raise valueerror("invalid gea method. valid methods include: %s" % (
                             ", ".join(methods)))
        return results


    ##############
    #getting data#
    ##############

    def get_coords(self, spp=0, individs=None):
        """
        returns the x,y coordinates of any or all individuals in a species

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        individs : iterable collection of ints, optional, default: None
            if provided, indicates the integer indices of the individuals
            whose coordinates should be returned.
            if None, all individuals' coordinates will be returned.

        returns
        -------
        np.array
            returns an nx2 numpy array, where n is the number of individuals
            whose coordinates were requested, and the 2 columns correspond to
            the x then y (i.e. lon then lat) coordinates.

        note: rows are sorted by index number of the individuals,
              regardless of input index order.
        """
        spp = self.comm[self._get_spp_num(spp)]
        if individs is not None:
            individs = np.sort(individs)
        coords = spp._get_coords(individs=individs)
        return coords


    def get_x(self, spp=0, individs=None):
        """
        returns the x coordinates of any or all individuals in a species

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        individs : iterable collection of ints, optional, default: None
            if provided, indicates the integer indices of the individuals
            whose coordinates should be returned.
            if None, all individuals' coordinates will be returned.

        returns
        -------
        np.array
            returns an nx1 numpy array, where n is the number of individuals
            whose coordinates were requested.

        note: rows are sorted by index number of the individuals,
            regardless of input index order.
        """
        spp = self.comm[self._get_spp_num(spp)]
        if individs is not None:
            individs = np.sort(individs)
        x = spp._get_x(individs=individs)
        return x


    def get_y(self, spp=0, individs=None):
        """
        returns the y coordinates of any or all individuals in a species

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        individs : iterable collection of ints, optional, default: None
            if provided, indicates the integer indices of the individuals
            whose coordinates should be returned.
            if None, all individuals' coordinates will be returned.

        returns
        -------
        np.array
            returns an nx1 numpy array, where n is the number of individuals
            whose coordinates were requested.

        note: rows are sorted by index number of the individuals,
            regardless of input index order.
        """
        spp = self.comm[self._get_spp_num(spp)]
        if individs is not None:
            individs = np.sort(individs)
        y = spp._get_y(individs=individs)
        return y


    def get_cells(self, spp=0, individs=None):
        """
        returns the cell coordinates of any or all individuals in a species

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        individs : iterable collection of ints, optional, default: None
            if provided, indicates the integer indices of the individuals
            whose cell coordinates should be returned.
            if None, all individuals' cell coordinates will be returned.

        returns
        -------
        np.array
            returns an nx2 numpy array, where n is the number of individuals
            whose coordinates were requested, and the 2 columns correspond to
            the i then j (i.e. lat then lon) array coordinates.

            note: ccolumns are in the reverse order of the columns
            of the array returned by `model.get_coords()`, to facilitate
            array-subsetting of the landscape layer raster arrays.

            note: rows are sorted by index number of the individuals,
            regardless of input index order.


        """
        spp = self.comm[self._get_spp_num(spp)]
        if individs is not None:
            individs = np.sort(individs)
        cells = spp._get_cells(individs=individs)
        return cells


    def get_random_individs(self, n, spp=0):
        """
        returns a size-n sample of index numbers
        for random individuals in a species.

        parameters
        ----------
        n : int
            the number of random individuals to sample

        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        returns
        -------
        np.array
            returns a length-n 1d numpy array containing the sorted index
            numbers of the sampled individuals.

        """
        spp = self.comm[self._get_spp_num(spp)]
        individs = np.sort(spp._get_random_individuals(n))
        return individs


    def get_e(self, spp=0, lyr_num=None, individs=None):
        """
        returns the current environmental values of any or all
        individuals in a species, for one or all landscape layers

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        individs : iterable collection of ints, optional, default: None
            if provided, indicates the integer indices of the individuals
            whose environmental values should be returned.
            if None, all individuals' environmental values  will be returned.

        lyr_num : int, default:None
            if provided, indicates the layer for which values should be
            returned. if None, values will be returned for all layers.

        returns
        -------
        np.array
            returns an nxl numpy array, where n is the number of individuals
            whose coordinates were requested, and l is the number of layers
            requested, ordered by increasing layer number.

            note: rows are sorted by index number of the individuals,
            regardless of input index order.
        """
        spp = self.comm[self._get_spp_num(spp)]
        if individs is not None:
            individs = np.sort(individs)
        e = spp._get_e(individs=individs, lyr_num=lyr_num)
        return e


    def get_z(self, spp=0, trt=None, individs=None):
        """
        returns the phenotypic values of any or all individuals in a species,
        for one or all traits

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        individs : iterable collection of ints, optional, default: None
            if provided, indicates the integer indices of the individuals
            whose phenotypes should be returned.
            if None, all individuals' phenotypes will be returned.

        trt : {int, str}, optional, default: 0
            a reference to the trait for which the gea should be run.
            can be either the trait's index number (i.e. its integer key in the
            genomicarchitecture's traits dict), or its name (as a character
            string).

        returns
        -------
        np.array
            returns an nxt numpy array, where n is the number of individuals
            whose coordinates were requested, and t is the number of traits
            requested.

            note: rows are sorted by index number of the individuals,
            regardless of input index order.
        """
        spp = self.comm[self._get_spp_num(spp)]
        assert (spp.gen_arch is not None and
                spp.gen_arch.traits is not None), ("species without traits "
                                                   "have no phenotypic values "
                                                   "to be calculated.")
        if trt is not None:
            trait_num = self._get_trt_num(spp, trt)
        else:
            trait_num = None
        if individs is not None:
            individs = np.sort(individs)
        z = spp._get_z(individs=individs, trait_num=trait_num)
        return z


    def get_fitness(self, spp=0, trt=None, individs=None):
        """
        returns the current fitness of any or all individuals in a species,
        for one or all traits

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        individs : iterable collection of ints, optional, default: None
            if provided, indicates the integer indices of the individuals
            whose fitness should be returned.
            if None, all individuals' fitness values will be returned.

        trt : {int, str}, optional, default: 0
            a reference to the trait for which the gea should be run.
            can be either the trait's index number (i.e. its integer key in the
            genomicarchitecture's traits dict), or its name (as a character
            string).

        returns
        -------
        np.array
            returns an nxt numpy array, where n is the number of individuals
            whose coordinates were requested, and t is the number of traits
            requested.

            note: rows are sorted by index number of the individuals,
            regardless of input index order.

        """
        spp = self.comm[self._get_spp_num(spp)]
        assert (spp.gen_arch is not None and
                spp.gen_arch.traits is not None), ("species without traits "
                                                   "have no fitness values "
                                                   "to be calculated.")
        trait_num = self._get_trt_num(spp, trt)
        if individs is not None:
            individs = np.sort(individs)
        z = spp._get_z(individs=individs, trait_num=trait_num)
        fit = spp._get_fit(individs=individs)
        return fit


    def get_genotypes(self, spp=0, loci=None, individs=None, biallelic=False):
        """
        returns the genotypes of any or all individuals in a species
        and for any or all loci

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        loci : iterable collection of ints, optional, default: None
            if provided, indicates the loci for which genotypes should be
            returned. if None, genotypes will be returned for all loci.

        individs : iterable collection of ints, optional, default: None
            if provided, indicates the integer indices of the only individuals
            for whom genotypes should be returned.
            if None, all individuals' genotypes will be returned.

        biallelic : bool, optional, default: False
            if False, results are returned as nxl numpy arrays, where n is the
            number of individuals requested and l is the number of loci
            requested.
            if True, results are returned as nxlx2 numpy arrays, where the 2
            dimensions on the third axis represent the two haploid genotypes of
            each diploid individual.

        returns
        -------
        np.array
            returns a numpy array containing the requested genotypes.
            the shape of the array depends on the number of individuals and
            loci and the return format requested.
            the first and second axes represent the individuals and loci,
            returned in sorted numerical order. the optional third axis
            represents ploidy, if biallelic genotypes are requested.
            (see explanation of the argument 'biallelic', above,
            for more detail.)

            note: all axes are sorted by index number of the object they
            represent (individuals, loci), regardless of input index order.
        """
        spp = self.comm[self._get_spp_num(spp)]
        assert spp.gen_arch is not None, ("species without genomes "
                                          "have no genotypes.")
        if individs is not None:
            individs = np.sort(individs)
        if loci is not None:
            loci = np.sort(loci)

        gts = spp._get_genotypes(loci=loci, individs=individs,
                                 biallelic=biallelic)
        return gts


    def get_tree_sequence(self, spp=0):
        """
        returns the tskit.treesequence object for the indicated species

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).
        """
        spp = self.comm[self._get_spp_num(spp)]
        if not spp.gen_arch.use_tskit:
            raise exception(("treesequence object cannot be generated "
                             "for a species whose pedigree is not being "
                             "tracked by tskit."))
        assert spp.gen_arch is not None, ("species without genomes "
                                          "have no tskit data structures "
                                          "to be returned.")
        # prep the tablecollection, then get its treesequence
        spp._sort_and_simplify_table_collection()
        ts = spp._tc.tree_sequence()
        return ts


    ###############
    #altering data#
    ###############

    def remove_individuals(self,
                           spp=0,
                           n=None,
                           n_left=None,
                           individs=None,
                          ):
        """
        remove individuals from a species.

        parameters
        ----------
        spp : {int, str}, optional, default: 0
            a reference to the species for which values should be returned.
            can be either the species' index number (i.e. its
            integer key in the community dict), or its name (as a character
            string).

        n: {int}, optional, default: None
            the number of random individuals to be removed.
            If None then `n_out` or `individs` must be provided.

        n_left: {int}, optional, default: None
            the number of individuals left in the species
            after a random selection of individuals have been removed.
            if None then `n` or `individs` must be provided.

        individs: {list, tuple, numpy.ndarray}, optional, default: None
            an iterable containing the integer ids (i.e., the species dict's keys)
            of the individuals to be removed. if None then `n` or `n_out` must be provided.

        returns
        -------
        None
            alters the specified species in place by removing all specified
            (if `individs` is not None) or randomly drawn (if `n` is not None)
            individuals.

        """
        # get the desired species
        spp = self.comm[self._get_spp_num(spp)]

        # call that species' corresponding method
        spp._remove_individuals(individs=individs,
                                n=n,
                                n_left=n_left,
                                verbose=True,
                               )


    def add_individuals(self,
                        n,
                        coords,
                        recip_spp=0,
                        source_spp=None,
                        source_msprime_params=None,
                        individs=None,
                       ):
        """
        add individuals to a recipient species object, either using a second
        species object as the source population or feeding a dict of parameters
        into msprime to simulate the source population.

        parameters
        ----------
        n: int
            the number of individuals to be added.

        coords: {tuple, list, numpy.ndarray}
            a tuple, list of numpy.ndarray indicating the x,y coordinates where the
            new individuals should be added to the recipient species.
            if `coords` contains 2 values then all added individuals will be placed
            at the x,y coordinate pair indicated by those 2 values.
            otherwise, `coords` may be of size `n` x 2, indicating the x,y
            coordinate pair at which each of the `n` individuals will be placed.

        recip_spp: {gnx.structs.species.species, int, str}, optional, default:
            None
            an int or str referring to the species in this model to which
            individuals should be added.

        source_spp: {gnx.structs.species.species, int, str}, optional, default:
            None
            a geonomics.structs.species.species object to act as the source
            population from which individuals will be taken
            (or an int or str referring to that species within this model).

        source_msprime_params: dict, optional, default: None
            a dict of keyword args to be fed to `gnx.sim_msprime_individs()`,
            to parameterize the msprime source population model from which
            new individuals will be drawn. these include:
            - recomb_rate         (required)
            - mut_rate            (required)
            - demography          (optional; default: None)
            - population_size     (optional; default: None)
            - ancestry_model      (optional; default: None)
            - random_seed         (optional; default: None)
            see `gnx.sim_msprime_individs` for details.

        individs: {tuple, list, numpy.ndarray}, optional, default: None
            a 1d numpy.ndarray or equivalent array-like containing the int indices
            of the individuals to be taken from `source_spp` and added to
            `recip_spp`. If None, Indivduals from `source_spp` are randomly
            chosen.

        returns
        -------
        None
            alters `recip_spp` in place by adding the provided
            (if `source_spp` is not None) or simulated (if `source_msprime_params`
            is not None) individuals.

        """
        warnings.warn(("this function is new and has been checked but not "
                       "extensively tested. user beware. please immediately report "
                       "any suspected bugs at "
                       "https://github.com/erthward/geonomics/issues "
                       "or by email to drew.terasaki.hart@gmail.com"))

        # get the recipient species
        if isinstance(recip_spp, int) or isinstance(recip_spp, str):
            recip_spp = self.comm[self._get_spp_num(recip_spp)]

        # and make sure it has already been burned in
        assert recip_spp.burned, ("individuals cannot be manually added to a "
                                  "species that hasn't been burned in yet.")

        # make sure either a species object is given or a set of valid msprime args
        assert ((source_spp is not None and source_msprime_params is None) or
                (source_spp is None and source_msprime_params is not None)), (
                        "source population must be provided either as "
                        "a gnx species object ('source_spp') "
                        "or as a set of valid args for gnx.run_msprime_sim() "
                        "('source_msprime_params'), but not both.")
        if source_spp is not None:
            assert (isinstance(source_spp, species) or
                    isinstance(source_spp, int) or
                    isinstnace(source_spp, str)), ("the "
                        "value given to 'source_spp' identifies a gnx "
                        "species object that will serve as the source "
                        "population, so must be either a gnx species "
                        "instance (i.e., an object of the class "
                        "gnx.structs.species.species) whose individuals will "
                        "be added to the recipient species, or a valid int or str "
                        "identifier of another species in this model whose"
                        "individuals will be added to the recipeint species.")
            # get the source_spp from this model, if necessary
            if isinstance(source_spp, int) or isinstance(source_spp, str):
                source_spp = self.comm[self._get_spp_num(source_spp)]

        # call the species method
        recip_spp._add_individuals(n=n,
                                   coords=coords,
                                   land=self.land,
                                   source_spp=source_spp,
                                   source_msprime_params=source_msprime_params,
                                   individs=individs,
                                  )


    ##############
    #writing data#
    ##############

    def write_gendata(self, filepath, spp=0, n=None, include_fixed_sites=True):
        """
        A convenience function, designed for interactive use, that writes
        genetic data for a sample of n Individuals to the given filepath

        Parameters
        ----------
        filepath : str
            The path and filename indicating the destination of the output file.
            Output format ('VCF', 'FASTA') will be determined by the
            filepath extension ('.vcf', .'fasta').

        spp : {int, str}, optional, default: 0
            A reference to the Species for which values should be returned.
            Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string).

        n : {int, None}, optional, default: None
            The number of individuals of randomly sample.
            If None, output file will contain the entire current population.

        include_fixed_sites : bool, optional, default: True
            A flag indicating whether or not to include fixed sites. Only valid
            for VCF files. Defaults to True.

        Notes
        -----
        This function is only designed for use with small models
        or in interactive mode. For simulations being run to generate
        datasets for research use, it is preferable to generate a
        parameters file containing a section for data section
        (`gnx.make_parameters_file(..., data=True)`), then parameterize
        data collection there.
        """
        # get and check extension
        extension  = filepath.split('.')[-1].lower()
        assert extension in ['vcf', 'fasta'], ('Must provide valid file '
                                               'extension. Valid extensions '
                                               'include ".vcf" and ".fasta".')
        # get the species
        spp = self.comm[self._get_spp_num(spp)]
        # get the index-keyed dict of the sample of individuals
        sample = _get_adhoc_sample(spp, n)
        # get the genotypes
        gts = spp._get_genotypes(individs=[*sample], as_dict=True)
        # format and write out the file, per the extension
        if extension == 'vcf':
            gen_data = _format_vcf(sample, gts, spp.gen_arch,
                                   include_fixed_sites = include_fixed_sites)
        elif extension == 'fasta':
            gen_data = _format_fasta(sample, gts)

        # write the data out
        _write_file(filepath, gen_data)


    def write_geodata(self, filepath, spp=0, n=None):
        """
        A convenience function, designed for interactive use, that writes
        tabular geodata for a sample of n Individuals to the given filepath

        Parameters
        ----------
        filepath : str
            The path and filename indicating the destination of the output file.
            Output format ('CSV', 'Shapefile', 'GeoJSON') will be determined
            by the filepath extension ('.csv', '.shp', '.json').

        spp : {int, str}, optional, default: 0
            A reference to the Species for which values should be returned.
            Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string).

        n : {int, None}, optional, default: None
            The number of individuals of randomly sample.
            If None, output file will contain the entire current population.

        Notes
        -----
        This function is only designed for use with small models
        or in interactive mode. For simulations being run to generate
        datasets for research use, it is preferable to generate a
        parameters file containing a section for data section
        (`gnx.make_parameters_file(..., data=True)`), then parameterize
        data collection there.
        """
        # get and check extension
        extension  = filepath.split('.')[-1].lower()
        assert extension in ['csv', 'shp',
                             'json'], ('Must provide valid file '
                                       'extension. Valid extensions '
                                       'include ".csv", ".shp", and ".json".')
        # get the species
        spp = self.comm[self._get_spp_num(spp)]
        # get the sample of sorted individuals' index values
        sample = _get_adhoc_sample(spp, n)
        # format and write out the file, per the extension
        if extension == 'csv':
            _write_csv(filepath=filepath, individuals=sample)
        elif extension == 'shp':
            _write_shapefile(filepath=filepath, individuals=sample)
        elif extension == 'json':
            _write_geojson(filepath=filepath, individuals=sample)


    def write_tskit_table_collection(self, file_basename, spp=0, sep=','):
        """
        A convenience function, designed for interactive use, that writes the
        tskit.TableCollection containing the spatial pedigree into CSV tables

        Parameters
        ----------
        file_basename : str
            The path and file basename to use for the resulting table files.
            Each file will be appended with a capitalized suffix indicating
            the table it contains (NODES, EDGES, SITES, MUTATIONS, INDIVIDUALS).

        spp : {int, str}, optional, default: 0
            A reference to the Species for which values should be returned.
            Can be either the Species' index number (i.e. its
            integer key in the Community dict), or its name (as a character
            string).

        sep : str, optional, default: ','
            Character separator to use in the output CSV

        Notes
        -----
        This function is only designed for use with small models
        or in interactive mode. For simulations being run to generate
        datasets for research use, it is preferable to generate a
        parameters file containing a section for data section
        (`gnx.make_parameters_file(..., data=True)`), then parameterize
        data collection there.
        """
        # get the species
        spp = self.comm[self._get_spp_num(spp)]
        # list of the needed tables
        needed = ['nodes', 'edges', 'sites', 'mutations', 'individuals']
        # get the tskit.TableCollection as a dict of dicts, keeping only the
        # needed tables, and converting them in pandas.DataFrames
        tables = {k: _tskit_table_to_pandas(v) for (k,
                                v) in spp._tc.asdict().items() if k in needed}
        # write out each table
        for k, table in tables.items():
            filepath = file_basename + '_' + k.upper() + '.csv'
            table.to_csv(filepath, index=False, sep=sep)


