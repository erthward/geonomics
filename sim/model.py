#!/usr/bin/python
# model.py

'''
##########################################

Module name:          sim.model


Module contains:
                      - classes and functions to facilitate building a
                        Geonomics model and running multiple iterations
                        of it, including across multiple cores


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           07-06-18
Documentation:        URL


##########################################
'''

#geonomics imports
from structs import landscape, community, genome
from sim import burnin, data, stats
from utils import _str_repr_ as _sr_

#other imports
import numpy as np
import numpy.random as r
import random
from copy import deepcopy
import sys, os, traceback


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

#the Model class, which will contain the main model objects, build the
#burn-in and main functions, allow customization of the objects, manage and run
#iterations, and produce and save required stats and data
class Model:

    #######################
    ### SPECIAL METHODS ###
    #######################

    def __init__(self, name, params, verbose=False):

        #get the full input params dict
        self.params = deepcopy(params)
        #get the model params (to make the following code shorter)
        m_params = self.params.model

        #set the model name (which will come from the params filename)
        self.name = name

        #set the _verbose flag
        self._verbose = verbose
        #and set the private __term_width__  and __tab_len__ attributes
        self._set_term_width()
        self._set_tab_len()

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
        self.land = self._make_landscape()
        self.comm = self._make_community()

        #and make the self._never_run attribute False,
        #to prevent the land and community from being reset
        #on the very first iteration
        self._never_run = True

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
        #unless any population has a genomic architecture, indicating that its
        #genomes should be reassigned after burn-in; in that case it will be 
        #reset to False as soon as the genomes are reassigned
        self.reassign_genomes = None

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
        name_str = name_str % (_sr_._get_spacing(name_str), self.name)
        #get strings for the numbers of its, burn-in timesteps, main timesteps
        n_its_str = "Number of iterations:%s%i"
        n_its_str = n_its_str % (_sr_._get_spacing(n_its_str), self.n_its)
        burn_t_str= "Number of burn-in timesteps (minimum):%s%i"
        burn_t_str = burn_t_str % (_sr_._get_spacing(burn_t_str), self.burn_T)
        main_t_str = "Number of main timesteps:%s%i"
        main_t_str = main_t_str % (_sr_._get_spacing(main_t_str), self.T)
        #get strings for land and comm
        comm_str = "Populations:"
        pop_strs = ["%s %i: '%s'" % (_sr_._get_spacing(""),
            i, pop.name) for i, pop in self.comm.items()]
        pop_strs[0] = comm_str + pop_strs[0][len(comm_str):]
        comm_str = '\n'.join(pop_strs)
        land_str = "Layers:"
        lyr_strs = ["%s %i: '%s'" % (_sr_._get_spacing(""),
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
        stats_str = stats_str % (_sr_._get_spacing(stats_str),
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
        geo_data_str = geo_data_str % (_sr_._get_spacing(geo_data_str),
                                        ', '.join(geo_data_strs))
        gen_data_str = "Gen-data collected:%s{%s}"
        gen_data_str = gen_data_str % (_sr_._get_spacing(gen_data_str),
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


    #method for getting a population's number, given a name (str) or number
    def _get_pop_num(self, pop_id):
        if isinstance(pop_id, int):
            return(pop_id)
        elif isinstance(pop_id, str):
            #get nums for all pops with  matching names
            pop_nums = [k for k, pop in self.comm.items(
                                                    ) if pop.name == pop_id]
            assert len(pop_nums) == 1, ("Expected to find a single Population "
                "with a name matching the name provided (%s). Instead "
                "found %i.") % (pop_id, len(nums))
            pop_num = pop_nums[0]
            return pop_num
        else:
            raise ValueError(("The Population identifier must be either a str "
                "(indicating the Population's name) or an int (indicating "
                "its key in the Community dict). Instead, a %s was "
                "provided.") % str(type(pop_id)))

    #method for getting a Trait's number, give a name (str) or number
    def _get_trt_num(self, pop, trt_id):
        if isinstance(trt_id, int) or trt_id is None:
            return(trt_id)
        elif isinstance(trt, str):
            #get nums for all traits with matching names
            trt_nums = [k for k, trt in pop.gen_arch.traits.items(
                                                    ) if trt.name == trt_id]
            assert len(trt_nums) == 1, ("Expected to find a single Trait in "
                "Population '%s' with a name matching the name provided (%s). "
                "Instead found %i.") % (pop.name, trt_id, len(trt_nums))
            trt_num = trt_nums[0]
            return trt_num
        else:
            raise ValueError(("The Trait identifier must be either a str "
                "(indicating the Trait's name) or an int (indicating "
                "its key in the Population.gen_arch.traits dict). Instead, "
                "a %s was provided.") % str(type(trt_id)))


    #private method for determining the width of
    #the terminal on the current system
    def _set_term_width(self):
        self.__term_width__ = int(os.popen('stty size', 'r').read().split()[1])

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
                    [pop.gen_arch is not None for pop in self.comm.values()])

    #method to set seed (will be run when Model object is first created, if
    #called for in params)
    def _set_seeds(self):
        random.seed(self.seed)
        r.seed(self.seed)

    #method to wrap around landscape._make_landscape
    def _make_landscape(self):
        land = landscape._make_landscape(self.params)
        return(land)

    #a method to reset the land as necessary
    #(recopying or regeneratiing as necessary)
    def _reset_landscape(self, rand_landscape=True):
        #deepcopy the original land if necessary 
        #(else randomly generate new land)
        if not rand_landscape:
            #verbose output
            if self._verbose:
                print(('Copying the original Lanscape '
                       'for iteration %i...\n\n') % self.it)
            #deepcopy the original land
            self.land = deepcopy(self.orig_land)
            #and reset the land._changer.changes, if needed (so that they
            #point to the current land, not previous objects with possibly
            #updated attribute values)
            if self.land._changer is not None:
                #verbose output
                if self._verbose:
                    print(('Resetting the _LandscapeChanger.changes '
                                                        'object...\n\n'))
                self.land._changer._set_changes(self.land)
        else:
            #verbose output
            if self._verbose:
                print(('Creating new Landscape for '
                                    'iteration %i...\n\n') % self.it)
            self.land = self._make_landscape()


    #method to wrap around community._make_community
    def _make_community(self):
        comm = community._make_community(self.land, self.params, burn = True)
        return(comm)


    #a method to reset the community (recopying or regenerating as necessary)
    def _reset_community(self, rand_comm=True):
        if not rand_comm:
            #verbose output
            if self._verbose:
                print(('Copying the original community for '
                                    'iteration %i...\n\n') % self.it)
            self.comm = deepcopy(self.orig_comm)
            #and reset the pop._changer.changes objects for each pop,
            #if needed (so that they point to the current populations,
            #not previous ones with updated attribute values)
            for pop in self.comm.values():
                if pop._changer is not None:
                    #verbose output
                    if self._verbose:
                        print(('Resetting the pop._changer.changes '
                            'object for population " %s"...\n\n') % pop.name)
                    pop._changer._set_changes(pop)
        else:
            #verbose ouput
            if self._verbose:
                print(('Creating new community for '
                                    'iteration %i...\n\n') % self.it)
            self.comm = self._make_community()


    #method to make a data._DataCollector object for this model
    def _make_data_collector(self):
        data_collector = data._DataCollector(self.name, self.params)
        return data_collector


    #method to reset the self._data_collector attribute 
    #(a data._DataCollector object)
    def _reset_data_collector(self):
        self._data_collector = self._make_data_collector()


    #method to make a stats._StatsCollector object for this model
    def _make_stats_collector(self):
        stats_collector = stats._StatsCollector(self.name, self.params)
        return stats_collector


    #method to reset the self._stats_collector attribute 
    #(a stats._StatsCollector object)
    def _reset_stats_collector(self):
        self._stats_collector = self._make_stats_collector()


    #method to reset all the model's objects
    #(land, community, and associated) and attributes
    def _reset(self, rand_landscape=None, rand_comm=None, repeat_burn=None):
        #default to the self.rand_<_> attributes, if not otherwise provided
        if rand_landscape is None:
            rand_landscape = self.rand_landscape
        if rand_comm is None:
            rand_comm = self.rand_comm
        if repeat_burn is None:
            repeat_burn = self.repeat_burn

        #deepcopy the original landscape and comm if necessary
        #(if told not to randomize landscape and not at the beginning
        #of the initial iteration), else randomly generate new
        if not self._never_run:
            self._reset_landscape(rand_landscape)
            self._reset_community(rand_comm)
        else:
            self._never_run = False

        #reset the self.burn_t and self.t attributes to -1
        self._reset_t()
        if repeat_burn:
            self._reset_burn_t()

        #reset the community and population t attributes
        self.comm._reset_t()
        for pop in self.comm.values():
            pop._reset_t()

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
                print('Creating the burn-in function queue...\n\n')
            self.burn_fn_queue = self._make_fn_queue(burn = True)
        #verbose output
        if self._verbose:
            print('Creating the main function queue...\n\n')
        self.main_fn_queue = self._make_fn_queue(burn = False)


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
            queue.append(self.comm._set_t)
            for pop in self.comm.values():
                queue.append(pop._set_t)

        #append the set_age_stage methods to the queue
        for pop in self.comm.values():
            queue.append(pop._set_age_stage)
        #append the set_Nt methods
        for pop in self.comm.values():
            queue.append(pop._set_Nt)
        #append the do_movement_methods, if pop._move
        for pop in self.comm.values():
            if pop._move:
                queue.append(pop._do_movement)
        #append the do_pop_dynamics methods
        #FIXME: Consider whether the order of these needs to be specified, or
        #randomized, should people want to eventually simulate
        #multiple, interacting populations
        for pop in self.comm.values():
            queue.append(pop._do_pop_dynamics)

        #add the Changer.make_change, data._DataCollector._write_data, and 
        #stats._StatsCollector._write_stats methods, if this is not the burn-in
        #and if pop and/or land have Changer objects, or if the model has 
        #_DataCollector or _StatsCollector objects (in the self._data_collector 
        #and self._stats_collector attributes)
        if not burn:
            #add land._make_change method
            if self.land._changer is not None:
                queue.append(lambda: self.land._make_change(self.t))
            #add pop._make_change methods
            for pop in self.comm.values():
                if pop._changer is not None:
                    queue.append(pop._make_change)
            #add self.write_data method
            if self._data_collector is not None:
                queue.append(self.write_data)
            #add self.calc_stats method
            if self._stats_collector is not None:
                queue.append(self.calc_stats)

            #TODO depending how I integrate the Stats module, 
            #add stats functions to this queue too



        #add the burn-in function if need be
        if burn:
            #for pop in self.comm.values():
                #queue.append(lambda: pop.check_burned(self.burn_T))
            queue.append(lambda: self.comm._check_burned(burn_T = self.burn_T))
        return(queue)


    #method to generate timestep_verbose_output
    def _print_timestep_info(self, mode):
        verbose_msg = '%s:\t%i:%i\n' % (mode, self.it,
                                [self.burn_t if mode == 'burn' else self.t][0])
        pops_submsgs = ''.join(['\tPOP: %s%sN=%i\t(births=%i\tdeaths=%i)\n' %
                        (pop.name,
                        ' ' * (30 - len(pop.name)),
                        pop.Nt[:].pop(),
                        pop.n_births[:].pop(),
                        pop.n_deaths[:].pop()) for pop in self.comm.values()])
        verbose_msg = verbose_msg + pops_submsgs
        print(verbose_msg)
        print('\t' + '.' * (self.__term_width__ - self.__tab_len__))


    #method to run the burn-in or main function
    #queue (as indicated by mode argument)
    def _do_timestep(self, mode):
        #do a burn-in step
        if mode == 'burn':
            for fn in self.burn_fn_queue:
                if True not in [pop.extinct for pop in self.comm.values()]:
                    fn()
                else:
                    break
            if self._verbose:
                self._print_timestep_info(mode)
            #if the burn-in is complete, reassign the genomes if needed
            #and then set self.comm.burned = True
            if np.all([pop.burned for pop in self.comm.values()]):
                if self.reassign_genomes:
                    for pop in self.comm.values():
                        if pop.gen_arch is not None:
                            #verbose output
                            if self._verbose:
                                print(('Assigning genomes for '
                                    'population "%s"...\n\n') % pop.name)
                            genome._set_genomes(pop, self.burn_T, self.T)
                    #and then set the reassign_genomes attribute to False, so
                    #that they won'r get reassigned again during this iteration
                    self.reassign_genomes = False
                #mark the community as burned in
                self.comm.burned = True
                #verbose output
                if self._verbose:
                    print('Burn-in complete.\n\n')
        #or do a main step
        elif mode == 'main':
            for fn in self.main_fn_queue:
                if  True not in [pop.extinct for pop in self.comm.values()]:
                    fn()
                else:
                    break
            if self._verbose:
                self._print_timestep_info(mode)

        #then check if any populations are extinct and
        #return the correpsonding boolean
        extinct = np.any([pop.extinct for pop in self.comm.values()])
        #verbose output
        if extinct and self._verbose:
            print(('XXXX     Population %s went extinct. '
                'Iteration %i aborting.\n\n') % (' & '.join(
                ['"' + pop.name + '"' for pop in self.comm.values(
                                                ) if pop.extinct]), self.it))
        return(extinct)


    #method to set the model's next iteration
    def _set_next_iteration(self, core=None):
        #update the iteration numbers
        self._set_it()
        #verbose output
        if self._verbose:
            print('~' * self.__term_width__ + '\n\n')
            print('Setting up iteration %i...\n\n' % self.it)

        #reset the model (including the burn-in,
        #if self.repeat_burn or if this is the first iteration) 
        self._reset(rand_landscape = self.rand_landscape,
                         rand_comm = self.rand_comm,
                         repeat_burn = self.repeat_burn)


    #method for running the model's next iteration
    def _do_next_iteration(self):

        #set the next iteration to be run
        self._set_next_iteration()

        #loop over the burn-in timesteps running the burn-in
        #queue (if copied pop isn't already burned in)
        if (self.rand_comm
            or (not self.rand_comm and self.repeat_burn)
            or self.it == 0):
            #verbose output
            if self._verbose:
                print('Running burn-in, iteration %i...\n\n' % self.it)
            #until all populations have pop.burned == True
            while not np.all([pop.burned for pop in self.comm.values()]):
                #run a burn-in timestep
                extinct = self._do_timestep(mode = 'burn')
                #and end the iteration early if any population is extinct
                if extinct:
                    break

            #overwrite self.orig_comm with the burned-in community,
            #if necessary
            if not self.rand_comm and not self.repeat_burn:
                #verbose output
                if self._verbose:
                    print('Saving burned-in community...\n\n')
                self.orig_comm = deepcopy(self.comm)

        #verbose output
        if self._verbose:
            print('Running main model, iteration %i...\n\n' % self.it)
        #loop over the timesteps, running the run_main function repeatedly
        for t in range(self.T):
            #run a main timestep
            extinct = self._do_timestep('main')
            #and end the iteration early if any population is extinct 
            if extinct:
                break


        ################
        #public methods#
        ################

    #method to run the overall model
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
        subdirectory, containing further subdirectories for each Population.

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
        "GEONOMICS_params_<...>.py" file in the current working directory
        before `gnx.make_parameters_file` is called).

        >>> gnx.make_parameters_file()
        >>> mod = gnx.make_model()
        >>> mod.run(verbose = True)
        TODO: PUT TYPICAL MODEL OUTPUT HERE, EVEN THOUGH IT'S ONLY PRINTED?

        """

        #TODO: Add some handler for farming instances out to 
        #different nodes, if available?

        #set the self._verbose flag
        self._verbose = verbose

        #verbose output
        if self._verbose:
            print('\n\n' + '#' * self.__term_width__ + '\n\n')
            print('Running model "%s"...\n\n' % self.name)

        #loop over all the iterations
        while len(self.its) > 0:
            #do the next iteration
            try:
                self._do_next_iteration()
            except Exception as e:
                msg = ('XXXX\tAn error occurred during iteration %i, '
                    'timestep %i.\n \tPLEASE COPY AND REPORT THE FOLLOWING, '
                    'to help us debug Geonomics:\n')
                if self.comm.burned:
                    msg = msg % (self.it, self.t)
                else:
                    msg = msg % (self.it, self.burn_t)
                print(msg)
                print('<>' * int(self.__term_width__/2))
                print('Error message:\n\t%s\n\n' % e)
                traceback.print_exc(file=sys.stdout)
                print()
                print('<>' * int(self.__term_width__/2) + '\n\n')

        #verbose output
        if self._verbose:
            print('\n\nModel "%s" is complete.\n' % self.name)
            print('#' * self.__term_width__)

        #set the _verbose flag back to False
        self._verbose = False


    #method to run the model interactively from the command line; named 'walk'
    #to succinctly differentiate it from 'run'
    def walk(self, T=1, mode='main', verbose=True):
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
            If 'burn', it will run in burn-in mode. This entails the following:
                - Genomes will not yet be assigned, so genomic phenomena such
                  as crossing-over and natural selection will not occur;
                - Any scheduled landscape and demographic changes will not
                  occur;
                - Statiionarity statistics will be checked at the end of each
                  timestep, to determine if the Model has burned in (in which
                  case the model will stop running regardless of timestep,
                  genomes will be randomly assigned to all individuals, and
                  the Model will be deemed ready to run in main mode)
            If 'main', it will run in main mode (i.e. all functionalities will
            be used. (Note that for a Model to run in main mode it must first
            have been burned in.)
        verbose : bool, optional
            Whether or not to run the Model should provide written output.
            If True, formatted messages will be printed to STDOUT at each
            timestep. Defaults to True for `Model.walk`.

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

    TODO: DOUBLE CHECK THAT IT STILL WRITE ALL DATA AND STATISTICS THAT WOULD
BE EXPECTED WHEN RUN WITH Model.walk.

        Examples
        --------
        We can create and walk the default model (as long as there is no
        "GEONOMICS_params_<...>.py" file in the current working directory
        before `gnx.make_parameters_file` is called).

        >>> gnx.make_parameters_file()
        >>> mod = gnx.make_model()
        >>> #run the burn-in until it is complete
        >>> mod.walk(T = 1000, mode = 'burn')
    TODO: PUT TYPICAL MODEL OUTPUT HERE, EVEN THOUGH IT'S ONLY PRINTED?
        >>> #now run in main mode for 50 timesteps
        >>> mod.walk(T = 50, mode = 'main')
    TODO: PUT TYPICAL MODEL OUTPUT HERE, EVEN THOUGH IT'S ONLY PRINTED?

        """

        #throw an error and quit if mode == 'main' but the community is 
        #not burned in
        if mode == 'main' and not self.comm.burned:
            raise ValueError(("The Model.walk method cannot be run in 'main' "
                "mode if the Model's Community has not yet been burned in "
                "(i.e. if Model.comm.burned is False)."))

        #set the model's _verbose flag
        self._verbose = verbose

        #if verbose, add a space below the command line prompt, for readability
        if self._verbose:
            print('\n')

        #run the model for the stipulated number of timesteps
        for t in range(T):
            #exit if burn-in is complete
            if mode == 'burn' and self.comm.burned:
                #verbose output
                if self._verbose:
                    print('Burn-in complete.\n\n')
                break
            #reset mode, if mode is 'burn' and there is no mod.burn_fn_queue
            if mode == 'burn' and self.burn_fn_queue is None:
                #verbose output
                if self._verbose:
                    print(('No mod.burn_fn_queue was found. '
                          'Running mod.reset()...\n\n'))
                self._reset()
            extinct = self._do_timestep(mode = mode)
            #end the iteration early if any population is extinct
            if extinct:
                break
        #reset self._verbose to False
        self._verbose = False

    #method to use the self._data_collector object to sample and write data
    def write_data(self):
        self._data_collector._write_data(self.comm, self.it)

    #method to use the self._stats_collector object to sample and write data
    def calc_stats(self):
        self._stats_collector._calc_stats(self.comm, self.t, self.it)


    ##########
    #plotting#
    ##########

    #wrapper around Population._plot and Landscape._plot
    #TODO: allow pop_id to be a list of pop ids, to plot each pop in a
    #different color with a legend!
    def plot(self, pop=None, lyr=None, hide_land=False, individs=None,
            text=False, color='black', edge_color='face', text_color='black',
            colorbar=True, size=25, text_size=9, im_interp_method='nearest',
            land_cmap='terrain', pt_cmap=None, alpha=False,
             zoom_width=None, x=None, y=None, vmin=None, vmax=None):
        #get the lyr num
        lyr_num = self._get_lyr_num(layer)
        #if no pop_id, then call Landscape._plot
        if pop_id is None:
            self.land._plot(lyr_num=lyr_num, colorbar=colorbar, cmap=land_cmap,
                im_interp_method=im_interp_method, x=x, y=y,
                zoom_width=zoom_width, vmin=vmin, vmax=vmax)
        #or else plot the pop
        else:
            #get the pop
            pop = self.comm[self._get_pop_num(pop)]
            #feed args into pop._plot
            pop._plot(lyr_num=lyr_num, hide_land=hide_land,
                individs=individs, text=text, color=color, edge_color=edge_color, 
                text_color=text_color, colorbar=colorbar, size=size,
                text_size=text_size, im_interp_method=im_interp_method,
                land_cmap=land_cmap, pt_cmap=pt_cmap, alpha=alpha,
                zoom_width=zoom_width, x=x, y=y, vmin=vmin, vmax=vmax)
            #add pop name
            plt.suptitle(pop.name)

    #wrapper around Population._plot_density
    def plot_density(self, pop, normalize=False, individs=None,
            text=False, color='black', edge_color='face',
            text_color='black', size=25, text_size = 9,
            alpha=0.5, zoom_width=None, x=None, y=None):
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #feed args into pop._plot_density
        pop._plot_density(normalize=normalize, individs=individs, text=text,
            color=color, edge_color=edge_color, text_color=text_color,
            size=size, text_size=text_size, alpha=alpha,
            zoom_width=zoom_width, x=x, y=y)
        #add pop name
        plt.suptitle(pop.name)

    #wrapper around Population._plot_genotype
    def plot_genotype(self, pop, locus, lyr=None, by_dominance=False,
            individs=None, text=False, size=25, text_size = 9,
            edge_color='black', text_color='black', colorbar=True,
            im_interp_method='nearest', alpha=1, zoom_width=None, x=None,
            y=None):
        #get the lyr num
        lyr_num = self._get_lyr_num(layer)
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #feed args into pop._plot_genotype
        pop._plot_genotype(locus=locus, lyr_num=lyr_num, individs=individs, 
            text=text, size=size, text_size=text_size, edge_color=edge_color,
            text_color=text_color, colorbar=colorbar,
            im_interp_method=im_interp_method, alpha=alpha,
            by_dominance=by_dominance, zoom_width=zoom_width, x=x, y=y)
        #add pop name
        plt.suptitle(pop.name)

    #wrapper around Population._plot_phenotype
    #for a given trait
    def plot_phenotype(self, pop, trait, layer=None, individs=None, 
            text=False, size=25, text_size = 9, edge_color='black',
            text_color='black', colorbar=True, im_interp_method='nearest',
            alpha=1, zoom_width=None, x=None, y=None):
        #get the lyr num
        lyr_num = self._get_lyr_num(layer)       
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #return messages if population does not have genomes or traits
        if pop.gen_arch is None:
            print(("Model.plot_phenotype is not valid for populations "
                "without genomes.\n"))
            return
        elif pop.gen_arch.traits is None:
            print(("Model.plot_phenotype is not valid for populations "
                "without traits.\n"))
            return
        #get the trt_num
        trt_num = self._get_trt_num(pop, trait)
        #trt_num can't be None for plot_phenotype
        assert trt_num is not None, ("None is not a valid value for the "
            "'trait' arguemnt.")
        #feed args into pop._plot_phenotype
        pop._plot_phenotype(trait=trait, lyr_num=lyr_num, individs=individs, 
            text=text, size=size, text_size=text_size, edge_color=edge_color,
            text_color=text_color, colorbar=colorbar,
            im_interp_method=im_interp_method, alpha=alpha,
            zoom_width=zoom_width, x=x, y=y)
        #add pop name
        plt.suptitle(pop.name)

    #wrapper around Population._plot_fitness
    def plot_fitness(self, pop, trait=None, layer=None, individs=None, 
            text=False, size=25, text_size = 9, edge_color='black',
            text_color='black', fit_cmap='RdYlGn', colorbar=True,
            im_interp_method='nearest', alpha=1, zoom_width=None, x=None,
            y=None):
        #get the lyr num
        lyr_num = self._get_lyr_num(layer)       
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #return messages if population does not have genomes or traits
        if pop.gen_arch is None:
            print(("Model.plot_fitness is not valid for populations "
                "without genomes.\n"))
            return
        elif pop.gen_arch.traits is None:
            print(("Model.plot_fitness is not valid for populations "
                "without traits.\n"))
            return
        #get the trt_num, which CAN be None for plot_fitness
        trt_num = self._get_trt_num(pop, trait)
        #feed args into pop._plot_fitness
        pop._plot_fitness(trt_num=trt_num, lyr_num=lyr_num, individs=individs, 
            text=text, size=size, text_size=text_size, edge_color=edge_color,
            text_color=text_color, fit_cmap=fit_cmap, colorbar=colorbar,
            im_interp_method=im_interp_method, alpha=alpha,
            zoom_width=zoom_width, x=x, y=y)
        #add pop name
        plt.suptitle(pop.name)

    #wrapper around pop._plot_allele_frequencies
    def plot_allele_frequencies(self, pop):
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #call the fn
        pop._plot_allele_frequencies()

    #wrapper around pop._plot_hist_fitness
    def plot_hist_fitness(self, pop):
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #call the fn
        pop._plot_hist_fitness()

    #wrapper around pop._plot_direction_surface for _move_surf
    def plot_movement_surface(self, pop, style, x, y, zoom_width=8,
                            scale_fact=4.5, color='black', colorbar = True):
        self._plot_direction_surface(surf_type='move', pop=pop, style=style,
            x=x, y=y, zoom_width=zoom_width, scale_fact=scale_fact,
            color=color, colorbar=colorbar)

    #wrapper around pop._plot_direciton_surface for _disp_surf
    def plot_dispersal_surface(self, pop, style, x, y, zoom_width=8,
                            scale_fact=4.5, color='black', colorbar = True):
        self._plot_direction_surface(surf_type='move', pop=pop, style=style,
            x=x, y=y, zoom_width=zoom_width, scale_fact=scale_fact,
            color=color, colorbar=colorbar)

    #wrapper around pop._plot_direction_surface
    def _plot_direction_surface(self, surf_type, pop, style, x, y, 
        zoom_width=8, scale_fact=4.5, color='black', colorbar = True):

        """
        The 'style' argument can take the following values:
            'hist':
                Plot a classic histogram approximating the Von Mises
                distribution at the cell indicated by position x,y.
            'circle_hist':
                Plot a circular histogram approximating the Von Mises
                distribution at the cell indicated by position x,y;
                plot will be drawn inside the chosen cell on the
                _DirectionalitySurface raster.
            'circle_draws':
                Plot points on the unit circle, whose locations were
                drawn from the Von Mises distribution at the cell indicated
                by position x,y; plot will be drawn inside the chosen cell
                on the _DirectionalitySurface raster.
            'vect':
                Inside each cell of the _DirectionalitySurface raster, plot the mean
                direction vector of directions drawn from that cell's Von Mises
                distribution.

        """
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #call the fn
        pop._plot_direction_surface(surf_type=surf_type, style=style, x=x, y=y,
            zoom_width=zoom_width, scale_fact=scale_fact, color=color,
            colorbar=colorbar)

    #wrapper around pop._plot_demographic_pyramid
    def plot_demographic_pyramid(self, pop):
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #call the fn
        pop._plot_demographic_pyramid()
    
    #wrapper around pop._plot_pop_growth
    def plot_pop_growth(self):
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #call the fn
        pop._plot_pop_growth()

    #wrapper around pop._plot_demographic_changes
    def plot_demographic_changes(self):
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #call the fn
        pop._plot_demographic_changes()

    #wrapper around pop._plot_stat
    def plot_stat(self, stat):
        #get the pop
        pop = self.comm[self._get_pop_num(pop)]
        #call the fn
        pop._plot_stat(stat)

