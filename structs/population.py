#!/usr/bin/python 
# population.py

'''
##########################################

Module name:              population


Module contains:
                          - definition of the Population class
                          - function for creating a population of Individuals (incl. their genomes and associated data)
                          - associated functions


Author:                    Drew Ellison Hart
Email:                     drew.hart@berkeley.edu
Github:                    URL
Start date:                12-28-15
Documentation:             URL


##########################################
'''

#geonomics imports
from utils import viz, spatial as spt
from structs import genome, landscape, individual
from ops import movement, mating, selection, mutation, demography, change
from sim import burnin

#other imports
import numpy as np
from numpy import random as r
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from collections import Counter as C
from collections import OrderedDict as OD
from copy import deepcopy
from operator import itemgetter
from operator import attrgetter
import sys


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

class Population(OD):

    #######################
    ### SPECIAL METHODS ###
    #######################

    def __init__(self, name, inds, land, pop_params, genomic_architecture=None):

        #check the inds object is correct
        assert type(inds) in (OD, dict), "Argument inds must be of type dict or type collections.Ordered_Dict"
        assert list(set([i.__class__.__name__ for i in list(inds.values())])) == [ 'Individual'], "Argument inds must be a dictionary containing only instances of the individual.Individual class"

        #then add the individuals
        self.update(inds) # update with the input dict of all individuals, as instances of individual.Individual class

        #set other attributes
        self.name = str(name)
        self.land_dim = land.dim
        self.land = land
        self.it = None  # attribute to keep track of iteration number this pop is being used for
            # (optional; will be set by the iteration.py module if called)
        self.T = None #attribute to track total number of timesteps to be run (will be set by the Model object)
        self.t = -1  # attribute to keep of track of number of timesteps the population has evolved (starts at -1, 
                     #to indicate unrun, and so that first timestep will be set to 0 at beginning of timestep)
            # NOTE: This way, if model is stopped, changes are made, then it is run further,
            # this makes it easy to continue tracking from the beginning
        self.burned = False #will be switched to True when the population
                            #passes the burn-in tests
        self.extinct = False #will be switched to True if the population goes extinct

        self.start_N = len(self)
        self.N = None  # attribute to hold a landscape.Scape object of the current population density
        self.K = None  # attribute to hold an landscape.Scape object of the local carrying capacity (i.e. 'target' dynamic equilibrium population density)
        self.Nt = []  # list to record population size (appended each time
        # pop.increment_age_stage() is called)
        self.n_births = []  # tracker of number of births each time pop.do_mating is called
        self.n_deaths = []  # tracker of number of deaths each time demography.pop_dynamics is called
        #attributes for storing numpy arrays of all individuals' coordinates and cells, to avoid repeat compute time each turn
        self.coords = None
        self.cells = None

        #create empty attributes to hold spatial objects that will be created after the population is instantiated
        self.kd_tree = None
        #create empty attribute to hold the DensityGridStack
        self.dens_grids = None

        #create empty attribute for a spatial.MovementSurface
        #which may be created, depending on paramters
        self.move_surf = None

        #create an empty changer attribute, which will be reset if the parameters define changes for this pop
        self.changer = None

        #set the sex_ratio to 0.5 default (but this will be changed if a sex
        #ratio is provided in the params and it is not 1/1)
        self.sex_ratio = 0.5
        #grab all of the mating, mortality, and movement parameters as Population attributes
        for section in ['mating', 'mortality', 'movement']:
            if section in [*pop_params]:
                for att,val in pop_params[section].items():
                    #leave out the move_surf components, which will be handled separately
                    if not isinstance(val, dict):
                        #convert sex ratio to the probabilty of drawing a male
                        if att == 'sex_ratio':
                            val = val / (val + 1)
                        setattr(self, att, val)

        #set the GenomicArchitecture object
        self.gen_arch = genomic_architecture
        assert (self.gen_arch.__class__.__name__ == 'GenomicArchitecture'
                or self.gen_arch is None), ("The Population.gen_arch attribute "
                "must be an instance of the genome.GenomicArchitecture class "
                "or else None.")

        #set the selection attribute, to indicate whether or not
        #natural selection should be implemented for the population
        self.selection = (self.gen_arch is not None and 
            (self.gen_arch.mu_delet > 0 or self.gen_arch.traits is not None))

        #set the self.mutate attribute (a boolean indicating whether or not to enact mutation, which is True if gen_arch.mu_tot > 0
        self.mutate = (self.gen_arch is not None
                       and self.gen_arch.mu_tot is not None
                       and self.gen_arch.mu_tot > 0)

        #set the self.mut_log attribute, which dictates whether or not a mutation log should be written for this pop
        self.mut_log = None
        if 'gen_arch' in [*pop_params]:
            self.mut_log = pop_params.gen_arch.mut_log

        #create a couple getter functions, for use in getting all individs' coordinates and certain individs' coordinates
        self.__coord_attrgetter__ = attrgetter('x', 'y')
        self.__individ_coord_attrgetter__ = attrgetter('idx', 'x', 'y')


    #override the __deepcopy__ method, so that the population can be copy.deepcopy'd
    #(because otherwise this doesn't work for classes that inherit from #collections.OrderedDict)
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
    #NOTE: this is not really a great representation; the Python docs indicate that __repr__ should ideally
    #provide a representation could be used to recreate the object, but if that is not possible then it 
    #should at least provide a string of the form '<... some useful description ...>; I've attempted to do
    #the latter, inspired by a combo of what I've seen in a few other packages (pandas, netCDF4, osgeo) 
    #(though perhaps I should more closely consider how to handle the params in a more nuanced, precise way once I'm
    #done writing the codebase, because right now this leaves out params that have lists and dictionarires at
    #values, and also includes some parameters that are input in params.py but not all, and some params that
    #are created internally but not all
    def __str__(self):
        #get a string representation of the class
        type_str = str(type(self))
        #get a string representation of the first and last individuals
        inds_str = '\n%i Individuals:\n\t' % len(self)
        first_ind_str = OD(self.items()).__str__().split('), ')[0] + ')\n\t...\n\t'
        last_ind_str = ', '.join(self.items().__str__().split(', ')[-2:])
        inds_str = inds_str + first_ind_str + last_ind_str + '\n'
        #get a string representation of the first two and last two parameters
        params = sorted([str(k) + ': ' +str(v) for k,v in vars(self).items() if type(v) in (str, int, bool, float)], key = lambda x: x.lower())
        params_str = "Parameters:\n\t" + ',\n\t'.join(params[:2]) + ','
        params_str = params_str + '\n\t...\n\t'
        params_str = params_str + ',\n\t'.join(params[-2:])
        return '\n'.join([type_str, inds_str, params_str])

    def __repr__(self):
        repr_str = self.__str__()
        return repr_str


    #####################
    ### OTHER METHODS ###
    #####################

    # method to set self.K
    def set_K(self, K):  # NOTE: Requires a landscape.Scape instance
        self.K = K

    # method to set self.N
    def set_N(self, N):  # NOTE: Requires a landscape.Scape instance
        self.N = N

    # method to append current pop size to the pop.Nt list
    def set_Nt(self):
        self.Nt.append(len(self))

    #method to increment the self.t attribute (i.e. the timestep counter)
    def set_t(self):
        self.t += 1

    #method to reset the self.t attribute (i.e. the timestep counter)
    def reset_t(self):
        self.t = -1

    # method to increment all population's age by one (also adds current pop size to tracking array)
    def set_age_stage(self):
        # increment age of all individuals
        [ind.set_age_stage() for ind in self.values()];
    # method to move all individuals simultaneously, and sample their new habitats

    def do_movement(self):
        movement.move(self)
        self.set_habitat()
        self.set_coords_and_cells()

    # function for finding all the mating pairs in a population
    def find_mating_pairs(self):
        mating_pairs = mating.find_mates(self, self.land)
        return mating_pairs

    # function for executing mating for a population
    def do_mating(self, mating_pairs, burn=False):

        #draw the number of births for each pair, and append total births to self.n_births list
        n_births = mating.draw_n_births(len(mating_pairs), self.n_births_distr_lambda)
        total_births = sum(n_births)
        self.n_births.append(total_births)

        #create the offspring_ids
        next_offspring_key = next(reversed(self))+1
        offspring_keys = set(range(next_offspring_key, next_offspring_key + total_births))
        #copy the keys, for use in mutation.do_mutation()
        keys_list = [*offspring_keys]

        if not burn and self.gen_arch is not None:
            recomb_paths = self.gen_arch.recomb_paths.get_paths(total_births*2)
            new_genomes = mating.do_mating(self, mating_pairs, n_births, recomb_paths)

        for n_pair, pair in enumerate(mating_pairs):

            parent_centroid_x = (self[pair[0]].x + self[pair[1]].x)/2
            parent_centroid_y = (self[pair[0]].y + self[pair[1]].y)/2

            n_offspring = n_births[n_pair]
            n_gametes = 2 * n_offspring

            for n in range(n_offspring):

                #get the next offspring_key
                offspring_key = offspring_keys.pop()

                offspring_x, offspring_y = movement.disperse(self.land, parent_centroid_x, parent_centroid_y,
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
                        new_genome = np.array([0])
                    elif self.gen_arch is None:
                        new_genome = None
                    else:
                        new_genome = new_genomes[n_pair][n]
                else:
                    new_genome = None

                #create the new individual
                self[offspring_key] = individual.Individual(idx = offspring_key, age = age, new_genome = new_genome, x = offspring_x, y = offspring_y, sex = sex)

                #set new individual's phenotype (won't be set during burn-in, because no genomes assigned; won't be set if the population has not gen_arch)
                if (self.gen_arch is not None 
                    and self.gen_arch.traits is not None 
                    and not burn):
                    self[offspring_key].set_phenotype(self.gen_arch)

        # sample all individuals' habitat values, to initiate for offspring
        self.set_habitat()
        self.set_coords_and_cells()

        # do mutation if necessary
        if self.mutate:
            mutation.do_mutation(keys_list, self, log = self.mut_log)


    #method to do population dynamics
    def do_pop_dynamics(self):
        #implement selection, iff self.selection is True and the pop has
        #already been burned in
        with_selection = self.selection and self.burned
        burn = not self.burned
        #then carry out the pop_dynamics, with selection as set above, and save
        #result, which will be True iff pop has gone extinct
        extinct = demography.do_pop_dynamics(self.land, self, with_selection = with_selection, burn = burn)
        if extinct:
            #set self.extinct equal to True, so that the iteration will be ended
            self.extinct = extinct

    #method to make population changes
    def make_change(self):
        self.changer.make_change(self.t)

    #method to check if the population has gone extinct
    def check_extinct(self):
        return len(self) == 0

    #method to check if the population is burned in (i.e. if it passes the burn-in checks)
        #TODO: Add more and/or different tests? This is a pretty basic test so far,
        #based only on gross population size (not spationarity of spatial distribution)
    #def check_burned(self, burn_T):
    #    burnin_status = len(self.Nt) >= burn_T
    #    if burnin_status:
    #        adf_test = burnin.test_adf_threshold(self, burn_T)
    #        t_test = burnin.test_t_threshold(self, burn_T)
    #        burnin_status = adf_test and t_test
    #    self.burned = burnin_status

    #method to calculate population density
    def calc_density(self, normalize = False, as_landscape = False, set_N=False):

        '''
        Calculate an interpolated raster of local population density, using
        the spatial.DensityGridStack object stored in the 
        Population.dens_grids attribute.

        If normalize is True, the density raster will vary between 0 and 1 
        (after being normalized by its own maximum value). If false, it will
        vary between 0 and its maximum estimated density value. 
        '''
        #check validity of normalize argument
        assert type(normalize) is bool, ("The 'normalize' argument takes "
            "a boolean value.\n")
        #get population coordinates
        x = self.get_x_coords()
        y = self.get_y_coords()
        #calculate the density array
        dens = self.dens_grids.calc_density(x, y)
        #set min value to 0
        dens = np.clip(dens, a_min = 0, a_max = None)
        #normalize, if necessary 
        if normalize:
            # Use max_val to normalize the density raster to either 0 to its current max val or
            # 0 to 1, to make sure the interpolation didn't generate any values slightly outside this range
            norm_factor = dens.max() - dens.min()
            dens = (dens - dens.min()) / norm_factor
        #return as landscape, if necessary
        if as_landscape == True:
            dens = landscape.Scape(self.land_dim, dens)
        #set self.N attribute, if necessary
        if set_N:
            self.set_N(dens)
        #or else return the density array
        else:
            return dens

    #method to set the individuals' habitat attributes
    def set_habitat(self, individs = None):
        if individs is None:
            inds_to_set = self.values()
        else:
            ig = itemgetter(*individs)
            inds_to_set = ig(self)
            if type(inds_to_set) is individual.Individual:
                inds_to_set = (inds_to_set,)
        hab = [ind.set_habitat([scape.rast[int(ind.y), int(ind.x)] for scape in self.land.values()]) for ind in inds_to_set]

    #method to set the individuals' phenotype attributes 
    def set_phenotype(self):
        [ind.set_phenotype(self.gen_arch) for ind in self.values()];

    #method to set the individuals' fitness attributes
    def set_fitness(self, fit):
        [ind.set_fitness(f) for ind, f in zip(self.values(), fit)];
        

    #method to set population's coords and cells arrays
    def set_coords_and_cells(self):
        self.coords = self.get_coords()
        self.cells = np.int32(np.floor(self.coords))


    #method to set the population's kd_tree attribute (a spatial.KDTree)
    def set_kd_tree(self, leafsize = 100):
        self.kd_tree = spt.KDTree(coords = self.coords, leafsize = leafsize)


    #method to set the population's spatial.DensityGridStack attribute
    def set_dens_grids(self, widow_width = None):
        self.dens_grids = spt.DensityGridStack(land = self.land, window_width = self.dens_grid_window_width)


    # method to get individs' habitat values
    def get_habitat(self, scape_num=None, individs=None):
        if individs is None:
            if scape_num is None:
                habs = np.array([ind.habitat for ind in self.values()])
            else:
                habs = np.array([ind.habitat[scape_num] for ind in self.values()])
        else:
            ig = itemgetter(*individs)
            if scape_num is None:
                habs = {i:ind.habitat for i, ind in self.items()}
                habs = np.array(ig(habs))
            else:
                habs = {i:ind.habitat[scape_num] for i, ind in self.items()}
                habs = np.array(ig(habs))
        return habs

    def get_genotype(self, locus, biallelic=False, individs=None, 
                                                        by_dominance=True):

        if individs is None:
            individs = [*self]

        if biallelic:
            return {i: self[i].genome[locus, :] for i in [*self] if i in individs}

        else:
            if by_dominance == True:
                d = self.gen_arch.d[locus]
            else:
                d = 0.5
            return dict(zip(individs, self.gen_arch.dom_effects[h](
                np.array([ind.genome[locus,] for i, ind in self.items() if i in individs]))))

    #convenience method for getting a scalar attribute for some or all individs
    def get_scalar_attr(self, attr_name, individs=None):
        if individs is None:
            vals = np.array([getattr(ind, attr_name) for ind in self.values()])
        else:
            ig = itemgetter(*individs)
            vals = {i: getattr(ind, attr_name) for i, ind in self.items()}
            vals = np.array(ig(vals))
        return vals

    #convenience method for getting age of whole population
    def get_age(self, individs=None):
        ages = self.get_scalar_attr('age', individs = individs)
        return ages

    # convenience method for getting whole population's phenotype
    def get_phenotype(self, individs=None):
        zs = self.get_scalar_attr('phenotype', individs = individs)
        return zs

    #convenience method for getting whole population's fitnesses
    def get_fitness(self, individs = None):
        fits = self.get_scalar_attr('fitness', individs = individs)
        return fits

    def calc_fitness(self, trait_num = None, set_fit = True):
        fit = selection.calc_fitness(self, trait_num = trait_num)
        #set individuals' fitness attributes, if indicated
        if set_fit:
            self.set_fitness(fit)
        return fit

    def get_dom(self, locus):
        return {locus: self.gen_arch.h[locus]}

    def get_coords(self, individs = None, as_float = True):
        coords = list(map(self.__coord_attrgetter__, self.values()))
        if individs is not None: 
            ig = itemgetter(*individs)
            coords = ig(dict(zip([*self], coords)))

        if as_float == True:
            coords = np.float32(coords)
        else:
            coords = np.int32(np.floor(coords))

        #make sure it's at least 2d (in case a single individual is requested)

        coords = np.atleast_2d(coords)

        return coords


    def get_cells(self, individs = None):
        cells = self.get_coords(individs = individs, as_float = False)
        return cells


    def get_x_coords(self, individs=None):
        coords = self.get_coords(individs = individs)
        return coords[:,0]


    def get_y_coords(self, individs=None):
        coords = self.get_coords(individs = individs)
        return coords[:,1]

    #method to return an n-length list of random individs; return individuals, or indices, as indicated
    def get_random_individuals(self, n, return_format='index'):
        assert return_format in ['index', 'individual'], "Argument return_format can take only 'index' or 'individual' as values (defaults to 'index')."
        choices = choices = r.choice(list(range(len(self))), n)
        inds = np.array(list(self.keys()))[choices]
        if return_format=='individual':
            inds = [self[ind] for ind in inds]
        return(inds)

    #use the kd_tree attribute to find nearest neighbors either within the population, if within == True,
    #or between the population and the points provided, if within == False and points is not None
    def find_neighbors(self, dist, within=True, coords=None, k = 2):
        #NOTE: In lieu of a more sophisticated way of determining whether the kd_tree needs to be updated 
        #(i.e. if the population have undergone movement, mating, or mortality since the last time it was 
        #constructed), and in an effort to minimize the number of times it is constructed each time (since
        #it's not an inconsequential runtime, so telling the model to rebuild it after each move, birth, or
        #death step could be unncessarily costly), for now I am just telling the tree to be rebuilt each time 
        #the pop.find_neighbors() method is called!
        self.set_kd_tree()
        #if neighbors are to be found within the population, set coords to self.coords
        #(otherwise, the coords to find nearest neighbors with should have been provided)
        if within:
            coords = self.coords
        #query the tree
        dists, pairs = self.kd_tree.query(coords = coords, dist = dist, k = k)
        return(dists, pairs)

    #method for plotting the population (or a subset of its individuals, by ID) on top of a landscape (or landscape stack)
    def plot(self, scape_num=None, hide_land=False, individs=None, text=False, color='black',
            edge_color='face', text_color='black', colorbar=True, size=25, text_size=9, im_interp_method='nearest', 
            land_cmap='terrain', pt_cmap=None, alpha=False, zoom_width=None, x=None, y=None, vmin = None, vmax = None):
        #convert individs to a list (in case comes in as a numpy array)
        if individs is not None and not isinstance(individs, list):
            individs = list(individs)
        #get coords
        if individs is None:
            coords = self.coords
            if text:
                text = [*self]
        else:
            coords = self.get_coords(individs)
            if text:
                text = individs
        if not text:
            text = None
        #set the plt_lims
        plt_lims = viz.get_plt_lims(self.land, x, y, zoom_width)
        #plot the landscape(s)
        if hide_land:
            pass
        else:
            viz.plot_rasters(self.land, scape_num = scape_num, colorbar = colorbar, im_interp_method = im_interp_method, cmap = land_cmap, plt_lims = plt_lims)
        #and plot the individuals
        viz.plot_points(coords, scape_num = scape_num, color = color, edge_color = edge_color, text_color = text_color, size = size, text_size = text_size, alpha = alpha, text = text, plt_lims = plt_lims, pt_cmap = pt_cmap, vmin = vmin, vmax = vmax)


    #method for plotting the population on top of its estimated population-density raster
    def plot_density(self, normalize=False, individs=None, text=False, color='black', edge_color='face', 
            text_color='black', size=25, text_size = 9, alpha=0.5, zoom_width=None, x=None, y=None):
        assert type(normalize) is bool, ("The 'normalize' argument takes "
            "a boolean value.\n")
        dens = self.calc_density(normalize = normalize)
        plt_lims = viz.get_plt_lims(self.land, x, y, zoom_width)
        if normalize:
            viz.plot_rasters(dens, plt_lims = plt_lims)
        else:
            viz.plot_rasters(dens, plt_lims = plt_lims, vmax = dens.max())
        self.plot(hide_land=True, individs = individs, text = text, color=color, edge_color = edge_color, text_color = text_color,
               size=size, text_size = text_size, alpha=alpha, zoom_width = zoom_width, x = x, y = y)


    # method for plotting individuals colored by their genotype at a given locus
    def plot_genotype(self, locus, scape_num=None, individs=None, text=False, size=25, text_size = 9,
            edge_color='black', text_color='black', colorbar=True, im_interp_method='nearest',
            alpha=1, by_dominance=False, zoom_width=None, x=None, y=None):

        if by_dominance == True:
            genotypes = self.get_genotype(locus, by_dominance=True)
        else:
            genotypes = self.get_genotype(locus)

        if individs is not None:
            genotypes = {i:v for i,v in genotypes.items() if i in individs}

        #colors to match mpl.cmap 'terrain' palette extremes, but with hybrid a mix of the extremes rather than the yellow at the middle of the palette, for nicer viewing: blue = [0,0], light blue = [0,1], white = [1,1]
        colors = ['#3C22B4', '#80A6FF', '#FFFFFF']

        for n, genotype in enumerate([0.0, 0.5, 1.0]):
            genotype_individs = [i for i, g in genotypes.items() if np.atleast_1d(g)[0] == genotype]
            # plot if there are any individuals of this genotype
            if len(genotype_individs) >= 1:
                self.plot(scape_num = scape_num, individs = genotype_individs, text = text, color = colors[n], edge_color = edge_color, 
                    text_color = text_color, colorbar = colorbar, size = size, text_size = text_size, 
                    im_interp_method = im_interp_method, alpha = alpha, zoom_width = zoom_width, x = x, y = y,
                    vmin = 0, vmax = 1)


    # method for plotting individuals colored by their phenotypes for a given trait
    def plot_phenotype(self, trait, scape_num=None, individs=None, text=False, size=25, text_size = 9, 
            edge_color='black', text_color='black', colorbar=True, im_interp_method='nearest', 
            alpha=1, zoom_width=None, x=None, y=None):

        z = OD(zip([*self], self.get_phenotype()[:,trait]))
        if individs is not None:
            z = {i:v for i,v in z.items() if i in individs}

        self.plot(scape_num = scape_num, individs = individs, text = text, color = list(z.values()),
                pt_cmap = 'terrain', edge_color = edge_color, text_color = text_color, colorbar = colorbar, 
                size = size, text_size = text_size, im_interp_method = im_interp_method, alpha = alpha, 
                zoom_width = zoom_width, x = x, y = y, vmin = 0, vmax = 1)


    # method for plotting individuals colored by their overall fitnesses, or by their fitnesses for a single
    # trait (if trait is not None)
    def plot_fitness(self, scape_num=None, trait_num=None, individs=None, text=False, size=100, text_size = 9, 
            edge_color='black', text_color='black', fit_cmap = 'RdYlGn', colorbar=True, im_interp_method='nearest', 
            alpha=1, zoom_width=None, x=None, y=None):
        
        #return messages if population does not have genomes or traits
        if self.gen_arch is None:
            print(("Population.plot_fitness is not valid for populations "
                   "without genomes.\n"))
            return
        elif self.gen_arch.traits is None:
            print(("Population.plot_fitness is not valid for populations "
                   "without traits.\n"))
            return

        # get all individs' fitness values
        if trait_num is None:
            w = self.calc_fitness()
        else:
            w = self.calc_fitness(trait_num = trait_num)

        #filter out unwanted individs, if necessary
        w = OD(zip([*self], w))
        if individs is not None:
            w = {i:v for i,v in w.items() if i in individs}

        # calc minimum possible fitness (for phenotypes within 0 <= z <= 1, which in reality isn't a
        # constraint, but values lower than this will also be constrained to the minimum-value color for plotting)
        #NOTE: the np.atleast_2d(...).min() construct makes this work both for fixed and spatially varying phi
        if trait_num is None:
            min_fit = np.product([1 - np.atleast_2d(t.phi).min() for t in list(self.gen_arch.traits.values())])
        else:
            min_fit = 1 - np.atleast_2d(self.gen_arch.traits[trait_num].phi).min()

        #then get uneven cmap and cbar-maker (needs to be uneven to give color-resolution to values varying
        #between 1 and the minimum-fitness value assuming all phenotypes are constrained 0<=z<=1, but then also
        #allow gradually deepening reds for fitness values lower than that minimum value), using the min_fit val
        cmap, make_cbar_fn = viz.make_fitness_cmap_and_cbar_maker(min_val = min_fit, max_val = 1, cmap = fit_cmap, trait_num = trait_num)

        #plot the trait phenotype in larger circles first, if trait is not None
        if trait_num is not None:
            #plot the outer (phenotype) circles
            self.plot_phenotype(trait=trait_num, scape_num=scape_num, individs=individs, text=False,
                    size=size, text_size = text_size, edge_color=edge_color, text_color=text_color,
                    colorbar=colorbar, im_interp_method=im_interp_method, 
                    alpha=alpha, zoom_width=zoom_width, x=x, y=y)
            #make size smaller for the next layer of inner (fitness) circles
            size = round(0.2*size)

        self.plot(scape_num = scape_num, individs = individs, text = text, color = list(w.values()),
                pt_cmap = cmap, edge_color = edge_color, text_color = text_color, colorbar = colorbar, 
                size = size, text_size = text_size, im_interp_method = im_interp_method, alpha = alpha, zoom_width = zoom_width, x = x, y = y)

        #and make a colorbar for the fitness values 
        viz.make_fitness_cbar(make_cbar_fn, min_fit)

    #method to plot a population's allele frequencies
    def plot_allele_frequencies(self):
        if self.gen_arch is None:
            print(("Population.plot_allele_frequencies is not valid for "
                "populations without genomes.\n"))
        else:
            self.gen_arch.plot_allele_frequencies(self)


    def plot_hist_fitness(self):
        plt.hist(list(self.calc_fitness()))
        plt.xlabel('Fitness')
        plt.ylabel('Count')


    #method for plotting the movement surface (in various formats)
    def plot_movement_surface(self, style, x, y, zoom_width=8, scale_fact=4.5, color='black', colorbar = True):
        '''
        The 'style' argument can take the following values: 
            'hist': 
                Plot a classic histogram approximating the Von Mises 
                distribution at the cell indicated by position x,y.
            'circle_hist':
                Plot a circular histogram approximating the Von Mises
                distribution at the cell indicated by position x,y;
                plot will be drawn inside the chosen cell on the 
                MovementSurface raster.
            'circle_draws':
                Plot points on the unit circle, whose locations were 
                drawn from the Von Mises distribution at the cell indicated 
                by position x,y; plot will be drawn inside the chosen cell 
                on the MovementSurface raster.
            'vect':
                Inside each cell of the MovementSurface raster, plot the mean
                direction vector of directions drawn from that cell's Von Mises
                distribution.

        '''
        if self.move_surf is None:
            print('This Population appears to have no MovementSurface. Function not valid.')
            return
        elif style not in ['hist', 'circ_hist', 'vect', 'circ_draws']:
            print("The 'style' argument must take one of the following values: 'hist', 'circ_hist', 'vect', 'circ_draws'")
            return
        elif style == 'hist':
            plt.hist(r.choice(self.move_surf.surf[y,x,:], size = 10000, replace = True), bins=100, density=True, alpha=0.5)

        else:
            #display the movement-surface raster
            scape_num = self.move_surf.scape_num
            self.land[scape_num].plot(zoom_width = zoom_width, x = x, y = y)

            if style == 'circ_hist':
                v, a = np.histogram(r.choice(self.move_surf.surf[y,x,:], replace = True, size = 7500), bins=15)
                v = v / float(v.sum())
                a = [(a[n] + a[n + 1]) / 2 for n in range(len(a) - 1)]
                xs = [np.cos(a[n]) * 0.5 for n in range(len(a))]
                ys = [np.sin(a[n]) * 0.5 for n in range(len(a))]
                xs = np.array(xs) * v * scale_fact
                ys = np.array(ys) * v * scale_fact
                [plt.plot((x, (x + xs[n])), (y, (y + ys[n])), linewidth=2, color=color) for n in range(len(xs))]

            elif style == 'circ_draws':
                pts = [(np.cos(a), np.sin(a)) for a in r.choice(self.move_surf.surf[y,x,:], size = 1000, replace = True)]
                plt.scatter([pt[0] * 0.5 + x for pt in pts], [pt[1] * 0.5 + y for pt in pts], color='red', alpha=0.1, marker = '.')

            elif style == 'vect':
                def plot_one_cell(x, y):
                    # draw sample of angles from the Gaussian KDE representing the von mises mixture distribution (KDE)
                    samp = self.move_surf.surf[y,x,:]
                    # create lists of the x and y (i.e. cos and sin) components of each angle in the sample
                    x_vects = np.cos(samp)
                    y_vects = np.sin(samp)
                    # define the dx and dy distances used to the position the arrowhead
                    # (divide by sqrt(2)/2, to scale to the size of half of the diagonal of a cell)
                    dx = np.mean(x_vects) / np.sqrt(2)
                    dy = np.mean(y_vects) / np.sqrt(2)
                    # now plot the arrow
                    plt.arrow(x, y, dx, dy, alpha=0.75, color='black', head_width=0.24, head_length=0.32)

                # call the internally defined function as a nested list comprehension for all raster cells, which I believe should do its best to vectorize the whole operation
                [[plot_one_cell(j, i) for i in range(self.move_surf.surf.shape[0])] for j in range(self.move_surf.surf.shape[1])]


    # method for plotting a population pyramid
    def plot_demographic_pyramid(self):
        #make dict of female and male colors
        col_dict = {-1: 'cyan', 1: 'pink'}
        #create a figure
        fig = plt.figure()
        #variables to grab the max count
        max_count = 0
        #for each sex
        for sex_val in [*col_dict]:
            #get a counter
            counts = C([ind.age for ind in self.values() if ind.sex == int(sex_val < 0)])
            #grab the ages from it
            ages = [*counts]
            #and grab the counts from it (multiplying by -1 for females, to set one sex on either side of x=0, for the pyramid)
            counts = [sex_val*count for count in counts.values()]
            #update the max_count var
            max_count = max(max(counts), max_count)
            #then create the horizontal barplot
            plt.barh(ages, counts, color = col_dict[sex_val])
        #use max_count to set the x-limits and y-limits
        plt.xlim((-1*max_count-2, max_count+2))
        #set the axis labels
        plt.xlabel('Population (individuals)')
        plt.ylabel('Age (timesteps)')
        #then set the xlabels to positive numbers on either side
        locs, labels = plt.xticks()
        plt.xticks(locs, [str(int(loc)) for loc in np.abs(locs)])
        #add sex symbols as title
        plt.suptitle('\u2642%s\u2640' % ''.join([' ' for _ in range(20)]), size = 30)
        #show it
        plt.show()

    def plot_pop_growth(self):
        T = range(len(self.Nt))
        x0 = self.Nt[0] / self.K.sum()
        plt.plot(T, [demography.calc_logistic_soln(x0, self.R, t) * self.K.sum() for t in T], color='red')
        plt.plot(T, self.Nt, color='blue')
        plt.xlabel('t')
        plt.ylabel('N(t)')

    def plot_demographic_changes(self):
        if self.changer is None:
            print(("Population.plot_demographic_changes is not valid "
                "for populations with no PopulationChanger object.\n"))
 
        else:
            self.changer.plot_dem_changes(self)

    def write_pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)

######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#function that uses the params.pops[<pop_num>].init.['K_<>'] parameters 
#to make the initial carrying-capacity raster
def make_K(pop, land, K_scape_num, K_fact):
    K_rast = land[K_scape_num].rast * K_fact
    return K_rast


def make_population(land, name, pop_params, burn=False):
    #get pop's intializing params
    init_params = deepcopy(pop_params.init)

    #if this population should have genomes, create the genomic architecture
    if 'gen_arch' in pop_params.keys():
        g_params = pop_params.gen_arch
        #make genomic_architecture
        gen_arch = genome.make_genomic_architecture(pop_params = pop_params)
    else:
        gen_arch = None

    #make individs
    N = init_params.pop('N')
    #create an ordered dictionary to hold the individuals, and fill it up
    inds = OD()
    for idx in range(N):
        # use individual.create_individual to simulate individuals and add them to the population
        ind = individual.make_individual(idx = idx, offspring = False, dim = land.dim, genomic_architecture = gen_arch, burn = burn)
        inds[idx] = ind

    #create the population from those individuals
    pop = Population(name = name, inds = inds, land = land, pop_params = pop_params, genomic_architecture=gen_arch)

    #use the remaining init_params to set the carrying-capacity raster (K)
    pop.set_K(make_K(pop, land, **init_params))
    #set initial habitat values
    pop.set_habitat()
    #set initial coords and cells
    pop.set_coords_and_cells()
    #set the kd_tree
    pop.set_kd_tree()

    #set phenotypes, if the population has genomes
    if pop.gen_arch is not None and not burn:
        [ind.set_phenotype(pop.gen_arch) for ind in pop.values()]

    #make density_grid
    pop.set_dens_grids()

    #if movement and movement_surf
    if 'movement' in [*pop_params] and pop_params.movement.move:
        if 'move_surf' in pop_params.movement.keys():
            ms_params = deepcopy(pop_params.movement.move_surf)
            #grab the scape number for the scape that the movement surface is to be based on
            move_surf_scape_num = ms_params.pop('scape_num')
            #then grab the move_surf scape's raster
            move_rast = land[move_surf_scape_num].rast
            #make the movement surface and set it as the pop's move_surf attribute
            pop.move_surf = spt.MovementSurface(land[move_surf_scape_num], **ms_params)

    #if this population has changes parameterized, create a PopulationChanger object for it
    if 'change' in pop_params.keys():
        #grab the change params
        ch_params = pop_params.change
        #make PopulationChanger and set it to the pop's changer attribute
        if land.changer is not None:
            pop.changer = change.PopulationChanger(pop, ch_params, land = land)
        else:
            pop.changer = change.PopulationChanger(pop, ch_params, land = None)

    return pop


# function for reading in a pickled pop
def read_pickled_pop(filename):
    import cPickle
    with open(filename, 'rb') as f:
        pop = cPickle.load(f)
    return pop

