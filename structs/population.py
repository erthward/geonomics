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

#other imports
import numpy as np
from numpy import random as r
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import cKDTree
from collections import Counter as C
from collections import OrderedDict as OD
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

    def __init__(self, name, inds, land, pop_params, genomic_architecture):

        #check the inds object is correct
        assert type(inds) in (OD, dict), "ERROR: inds must be of type dict or type collections.Ordered_Dict"
        assert list(set([i.__class__.__name__ for i in list(inds.values())])) == [ 'Individual'], "ERROR: inds must be a dictionary containing only instances of the individual.Individual class"
        
        #then add the individuals
        self.update(inds) # update with the input dict of all individuals, as instances of individual.Individual class
        self.inds = self.values()

        #set other attributes
        self.name = name
        self.land_dim = land.dim
        self.it = None  # attribute to keep track of iteration number this pop is being used for
            # (optional; will be set by the iteration.py module if called)
        self.t = 0  # attribute to keep of track of number of timesteps the population has evolved
            # NOTE: This way, if model is stopped, changes are made, then it is run further,
            # this makes it easy to continue tracking from the beginning
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
        #create empty attribute to hold the Density_Grid_Stack
        self.dens_grids = None

        #create empty attribute for a spatial.Movement_Surface
        #which may be created, depending on paramters
        self.move_surf = None

        #create empty attributes for islands and related attributes,
        #which may be created, depending on paramters
        self.islands = None
        self.n_island_mask_scape = None

        #create an empty changer attribute, which will be reset if the parameters define changes for this pop
        self.changer = None

        #grab all of the mating, mortality, and movement parameters as Population attributes
        for section in ['mating', 'mortality', 'movement']:
            for att,val in pop_params[section].items():
                #leave out the move_surf and islands components, which will be handled separately
                if not isinstance(val, dict):
                    setattr(self, att, pop_params[section][att])

        #TODO: probably get rid of this, or move it elsewhere in the package?
        #set the heterozygote effects dictionary
        self.heterozygote_effects = {
            0: lambda g: [np.ceil(np.mean(g[i,])) for i in range(g.shape[0])],
            # relative fitness of homozygous 1 = relative fitness of heterozygote, i.e. 1 = 1-hs
            0.5: lambda g: [np.mean(g[i,]) for i in range(g.shape[0])],
            # relative fitness of heterozygote halfway between homozygotes, i.e 1-hs = 1-s/2
            1: lambda g: [np.floor(np.mean(g[i,])) for i in range(g.shape[0])]
            # relative fitness of homozygous 0 = relative fitness of heterozygote, i.e. 1-hs = 1-s
        }

        #set the Genomic_Architecture object
        self.gen_arch = genomic_architecture
        assert self.gen_arch.__class__.__name__ == 'Genomic_Architecture', "self.gen_arch must be an instance of the genome.Genomic_Architecture class"

        #create a couple getter functions, for use in getting all individs' coordinates and certain individs' coordinates
        self.__coord_attrgetter__ = attrgetter('x', 'y')
        self.__individ_coord_attrgetter__ = attrgetter('idx', 'x', 'y')

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
        type_str = str(type(self))
        inds_str = '\n%i Individuals:\n\t' % len(self)
        first_ind_str = OD(self.items()).__str__().split('), ')[0] + ')\n\t...\n\t'
        last_ind_str = ', '.join(self.items().__str__().split(', ')[-2:])
        inds_str = inds_str + first_ind_str + last_ind_str + '\n'
        params = sorted([str(k) + ': ' +str(v) for k,v in vars(self).items() if type(v) in (str, int, bool, float)], key = lambda x: x.lower()) 
        params_str = "Parameters:\n\t" + ',\n\t'.join(params) 
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
        #self.Nt.append(self.get_size())

    # method to increment all population's age by one (also adds current pop size to tracking array)
    def reset_age_stage(self, burn=False):

        # increment age of all individuals
        [ind.reset_age_stage() for ind in self.inds];

        # add 1 to pop.t
        if burn == False:  # as long as this is not during burn-in, pop.t will increment
            self.t += 1

    # method to move all individuals simultaneously, and sample their new habitats
    def do_movement(self, land):
        movement.move(self)
        self.set_habitat(land)
        self.set_coords_and_cells()

    # function for finding all the mating pairs in a population
    def find_mating_pairs(self, land):
        mating_pairs = mating.find_mates(self, land)
        return mating_pairs

    # function for executing mating for a population
    def do_mating(self, land, mating_pairs, burn=False):

        n_births = mating.draw_n_births(len(mating_pairs), self.n_births_lambda)
        total_births = sum(n_births)
        self.n_births.append(total_births)

        if burn == False:
            recomb_paths = self.gen_arch.recomb_paths.get_paths(total_births*2)
            new_genomes = mating.do_mating(self, mating_pairs, n_births, recomb_paths)
        
        for n_pair, pair in enumerate(mating_pairs):

            parent_centroid_x = (self[pair[0]].x + self[pair[1]].x)/2
            parent_centroid_y = (self[pair[0]].y + self[pair[1]].y)/2

            n_offspring = n_births[n_pair]
            n_gametes = 2 * n_offspring

            offspring_key = next(reversed(self)) + 1
            for n in range(n_offspring):

                offspring_x, offspring_y = movement.disperse(land, parent_centroid_x, parent_centroid_y,
                        self.dispersal_mu, self.dispersal_sigma)

                if self.sex:
                    offspring_sex = r.binomial(1, 0.5)

                age = 0

                if not burn:
                    if self.sex:
                        self[offspring_key] = individual.Individual(idx = offspring_key, new_genome = new_genomes[n_pair][n], x = offspring_x, y = offspring_y, sex = offspring_sex, age =age)
                    else:
                        self[offspring_key] = individual.Individual(idx = offspring_key, new_genome = new_genomes[n_pair][n], x = offspring_x, y = offspring_y, age = age)
               
                    #set new individual's phenotype (won't be set during burn-in, becuase no genomes assigned)
                    self[offspring_key].set_phenotype(self.gen_arch)


                elif burn:
                    if self.sex:
                        self[offspring_key] = individual.Individual(idx = offspring_key, new_genome = np.array([0]), x = offspring_x, y = offspring_y, sex = offspring_sex, age = age)
                    else:
                        self[offspring_key] = individual.Individual(idx = offspring_key, new_genome = np.array([0]), x = offspring_x, y = offspring_y, age = age)
                offspring_key += 1

        # sample all individuals' habitat values, to initiate for offspring
        self.set_habitat(land)
        self.set_coords_and_cells()

        print('\n\t%i individuals born' % (total_births))

    # method to carry out mutation
    def do_mutation(self, log=False):
        mutation.do_mutation(self)

    #method to make population changes
    def make_change(self):
        self.changer.make_change(self.t)

    #method to check if the population has gone extinct
    def check_extinct(self):
        if len(self) == 0:
            print('\n\nYOUR POPULATION WENT EXTINCT.\n\n')
            return 1
            # sys.exit()
        else:
            return 0

    def calc_density(self, normalize_by= None, min_0=True, max_1=False, max_val=None, as_landscape = False, set_N=False):

        '''
        Calculate an interpolated raster of local population density, using spatial.Density_Grid_Stack object
        in self.dens_grids.

        Valid values for normalize_by currently include 'pop_size' and 'none'. If normalize_by = 'pop_size', max_1 =
        True will cause the output density raster to vary between 0 and 1, rather than between 0 and the current
        max normalized density value. 
        '''

        x = self.get_x_coords()
        y = self.get_y_coords()

        dens = self.dens_grids.calc_density(x, y)

        if normalize_by is not None:

            # if max_1 == True, set max_val to dens.max(), such that the density raster output will be normalized to
            # its own max, and thus vary between 0 and 1; else set to 1, and the output raster will vary between 0 and the current max value
            if max_1 == True:
                max_val = dens.max()
            elif max_1 == False:
                max_val = 1

            # Use max_val to normalize the density raster to either 0 to its current max val or
            # 0 to 1, to make sure the interpolation didn't generate any values slightly outside this range
            norm_factor = max_val - dens.min()
            dens = (dens - dens.min()) / norm_factor

        if min_0:
            dens[dens < 0] = 0

        if max_val is not None:
            dens[dens > max_val] = max_val

        if as_landscape == True:
            dens = landscape.Scape(land.dim, dens)

        if set_N:
            self.set_N(dens)

        else:
            return dens


    def set_habitat(self, land, individs = None):
        if individs is None:
            inds_to_set = self.inds
        else:
            ig = itemgetter(*individs)
            inds_to_set = ig(self)
            if type(inds_to_set) is individual.Individual:
                inds_to_set = (inds_to_set,)
        hab = [ind.set_habitat([scape.rast[int(ind.y), int(ind.x)] for scape in land.scapes]) for ind in inds_to_set]

    
    def set_phenotype(self):
        [ind.set_phenotype(self.gen_arch) for ind in self.inds];


    def set_fitness(self):
        indices = [*self]  #get list of individs' indices
        fit = self.get_fitness()
        [self[i].set_fitness(fit[n]) for n, i in enumerate(indices)]
        

    #method to set population's coords and cells arrays
    def set_coords_and_cells(self):
        self.coords = self.get_coords()
        self.cells = np.int32(np.floor(self.coords))


    #method to set the population's spatial.KD_Tree attribute
    def set_kd_tree(self, leafsize = 100):
        self.kd_tree = spt.KD_Tree(coords = self.coords, leafsize = leafsize)


    #method to set the population's spatial.Density_Grid_Stack attribute
    def set_dens_grids(self, land, widow_width = None):
        self.dens_grids = spt.Density_Grid_Stack(land = land, window_width = self.dens_grid_window_width)


    # method to get individs' habitat values
    def get_habitat(self, scape_num=None, individs=None):
        if individs is None:
            if scape_num is None:
                habs = np.array([ind.habitat for ind in self.inds])
            else:
                habs = np.array([ind.habitat[scape_num] for ind in self.inds])
        else:
            ig = itemgetter(*individs)
            if scape_num is None:
                habs = {i:ind.habitat for i, ind in self.items()}
                habs = np.array(ig(habs))
            else:
                habs = {i:ind.habitat[scape_num] for i, ind in self.items()}
                habs = np.array(ig(habs))
        return habs

    def get_age(self, individs=None):
        if individs is None:
            ages = np.array([ind.age for ind in self.inds])
        else:
            ig = itemgetter(*individs)
            ages = {i: ind.age for i, ind in self.items()}
            ages = np.array(ig(ages))
        return ages

    def get_genotype(self, locus, return_format='mean', individs=None, by_dominance=False):

        if individs is None:
            individs = [*self]

        if return_format == 'biallelic':
            return {i: self[i].genome[locus, :] for i in [*self] if i in individs}

        elif return_format == 'mean':
            if by_dominance == True:
                h = self.gen_arch.h[locus]
            else:
                h = 0.5
            return dict(zip(individs, self.heterozygote_effects[h](
                np.array([ind.genome[locus,] for i, ind in self.items() if i in individs]))))

    # convenience method for getting whole population's phenotype
    def get_phenotype(self, individs=None):
        if individs is None:
            zs = np.array([ind.phenotype for ind in self.inds])
        else: 
            ig = itemgetter(*individs)
            zs = {i:ind.phenotype for i, ind in self.items()}
            zs = np.array(ig(zs))
        return zs

    def get_fitness(self, trait_num = None, set_fit = False):
        return selection.calc_fitness(self, trait_num = trait_num, set_fit = set_fit)

    def get_dom(self, locus):
        return {locus: self.gen_arch.h[locus]}

    def get_coords(self, individs = None, as_float = True):
        coords = list(map(self.__coord_attrgetter__, self.inds))
        if individs is not None: 
            ig = itemgetter(*individs)
            coords = ig(dict(zip([*self], coords)))

        if as_float == True:
            coords = np.float32(coords)
        else:
            coords = np.int32(np.floor(coords))

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
    def show(self, land, scape_num=None, hide_land=False, individs=None, text=False, color='black',
            edge_color='face', text_color='black', colorbar=True, size=25, text_size=9, im_interp_method='nearest', 
            land_cmap='terrain', pt_cmap=None, alpha=False, zoom_width=None, x=None, y=None, vmin = None, vmax = None):
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
        plt_lims = viz.get_plt_lims(land, x, y, zoom_width)
        #plot the landscape(s)
        if hide_land:
            pass
        else:
            viz.show_rasters(land, scape_num = scape_num, colorbar = colorbar, im_interp_method = im_interp_method, cmap = land_cmap, plt_lims = plt_lims)
        #and plot the individuals
        viz.show_points(coords, scape_num = scape_num, color = color, edge_color = edge_color, text_color = text_color, size = size, text_size = text_size, alpha = alpha, text = text, plt_lims = plt_lims, pt_cmap = pt_cmap, vmin = vmin, vmax = vmax)


    #method for plotting the population on top of its estimated population-density raster
    def show_density(self, land, normalize_by='pop_size', individs=None, text=False, max_1=False, color='black', edge_color='face', 
            text_color='black', size=25, text_size = 9, alpha=0.5, zoom_width=None, x=None, y=None):
        dens = self.calc_density(normalize_by=normalize_by, max_1=max_1)
        plt_lims = viz.get_plt_lims(land, x, y, zoom_width)
        viz.show_rasters(dens, plt_lims = plt_lims)
        self.show(land, hide_land=True, individs = individs, text = text, color=color, edge_color = edge_color, text_color = text_color, 
               size=size, text_size = text_size, alpha=alpha, zoom_width = zoom_width, x = x, y = y)


    # method for plotting individuals colored by their genotype at a given locus
    def show_genotype(self, locus, land, scape_num=None, individs=None, text=False, size=25, text_size = 9, 
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
                self.show(land, scape_num = scape_num, individs = genotype_individs, text = text, color = colors[n], edge_color = edge_color, 
                    text_color = text_color, colorbar = colorbar, size = size, text_size = text_size, 
                    im_interp_method = im_interp_method, alpha = alpha, zoom_width = zoom_width, x = x, y = y,
                    vmin = 0, vmax = 1)


    # method for plotting individuals colored by their phenotypes for a given trait
    def show_phenotype(self, trait, land, scape_num=None, individs=None, text=False, size=25, text_size = 9, 
            edge_color='black', text_color='black', colorbar=True, im_interp_method='nearest', 
            alpha=1, by_dominance=False, zoom_width=None, x=None, y=None):

        z = OD(zip([*self], self.get_phenotype()[:,trait]))
        if individs is not None:
            z = {i:v for i,v in z.items() if i in individs}

        self.show(land, scape_num = scape_num, individs = individs, text = text, color = list(z.values()),
                pt_cmap = 'terrain', edge_color = edge_color, text_color = text_color, colorbar = colorbar, 
                size = size, text_size = text_size, im_interp_method = im_interp_method, alpha = alpha, 
                zoom_width = zoom_width, x = x, y = y, vmin = 0, vmax = 1)


    # method for plotting individuals colored by their overall fitnesses, or by their fitnesses for a single
    # trait (if trait is not None)
    def show_fitness(self, land, scape_num=None, trait_num=None, individs=None, text=False, size=100, text_size = 9, 
            edge_color='black', text_color='black', fit_cmap = 'RdYlGn', colorbar=True, im_interp_method='nearest', 
            alpha=1, by_dominance=False, zoom_width=None, x=None, y=None):

        # get all individs' fitness values
        if trait_num is None:
            w = self.get_fitness()
        else:
            w = self.get_fitness(trait_num = trait_num)

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
        cmap, make_cbar = viz.make_fitness_cmap_and_cbar_maker(min_val = min_fit, max_val = 1, cmap = fit_cmap, trait_num = trait_num)
        
        #plot the trait phenotype in larger circles first, if trait is not None
        if trait_num is not None:
            #plot the outer (phenotype) circles
            self.show_phenotype(trait=trait_num, land=land, scape_num=scape_num, individs=individs, text=False,
                    size=size, text_size = text_size, edge_color=edge_color, text_color=text_color,
                    colorbar=colorbar, im_interp_method=im_interp_method, 
                    alpha=alpha, by_dominance=by_dominance, zoom_width=zoom_width, x=x, y=y)
            #make size smaller for the next layer of inner (fitness) circles
            size = round(0.2*size)

        self.show(land, scape_num = scape_num, individs = individs, text = text, color = list(w.values()),
                pt_cmap = cmap, edge_color = edge_color, text_color = text_color, colorbar = colorbar, 
                size = size, text_size = text_size, im_interp_method = im_interp_method, alpha = alpha, zoom_width = zoom_width, x = x, y = y)

        #and make a colorbar for the fitness values 
        viz.make_fitness_cbar(make_cbar, min_fit)


    def show_hist_fitness(self):
        plt.hist(list(selection.get_fitness(self).values()))
 

    # method for plotting a population pyramid
    # NOTE: NEED TO FIX THIS SO THAT EACH HALF PLOTS ON OPPOSITE SIDES OF THE Y-AXIS
    def show_pyramid(self):
        plt.hist([ind.age for ind in list(self.inds) if ind.sex == 0], orientation='horizontal',
                 color='pink',
                 alpha=0.6)
        plt.hist([ind.age for ind in list(self.inds) if ind.sex == 1], orientation='horizontal',
                 color='skyblue',
                 alpha=0.6)


    def show_pop_growth(self):
        T = range(len(self.Nt))
        x0 = self.Nt[0] / self.K.sum()
        plt.plot(T, [demography.logistic_soln(x0, self.R, t) * self.K.sum() for t in T], color='red')
        plt.plot(T, self.Nt, color='blue')
        plt.xlabel('t')
        plt.ylabel('N(t)')


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


def make_population(land, pop_params, burn=False):
    #get pop's intializing params
    init_params = pop_params.init

    #if this population should have genomes, create the genomic architecture
    if 'genome' in pop_params.keys():
        g_params = pop_params.genome
        #make genomic_architecture
        gen_arch = genome.make_genomic_architecture(g_params = g_params)

    #make individs
    N = init_params.pop('N')
    #create an ordered dictionary to hold the individuals, and fill it up
    inds = OD()
    for idx in range(N):
        # use individual.create_individual to simulate individuals and add them to the population
        ind = individual.make_individual(idx = idx, offspring = False, dim = land.dim, genomic_architecture = gen_arch, burn = burn)
        inds[idx] = ind
    
    #create the population from those individuals
    name = init_params.pop('name')
    pop = Population(name = name, inds = inds, land = land, pop_params = pop_params, genomic_architecture=gen_arch)
  
    #use the remaining init_params to set the carrying-capacity raster (K)
    pop.set_K(make_K(pop, land, **init_params))
    #set initial habitat values
    pop.set_habitat(land)
    #set initial coords and cells
    pop.set_coords_and_cells()
    #set the kd_tree
    pop.set_kd_tree()

    #set phenotypes, if the population has genomes
    if pop.gen_arch is not None and not burn:
        [ind.set_phenotype(pop.gen_arch) for ind in pop.inds]

    #make density_grid
    pop.set_dens_grids(land = land)

    #if movement and movement_surf
    if pop_params.movement.move:
        if 'move_surf' in pop_params.movement.keys():
            ms_params = pop_params.movement.move_surf
            #TODO: Once I have a meta-function that takes arguments about whether or not a movement_surface is
            #desired for each pop, and the same for islands, and then generates a template params file, then I
            #will get rid of the 'make' keys here and in the islands section below
            make_ms = ms_params.pop('make')
            if make_ms:
                #grab the scape number for the scape that the movement surface is to be based on
                move_surf_scape_num = ms_params.pop('scape_num')
                #then grab the move_surf scape's raster
                move_rast = land[move_surf_scape_num].rast

###ISLANDS############################################################################################################
                #TODO: EITHER GET RID OF, OR DEEPLY RETHINK, THE ISLANDS THING
                #make the islands layer if required
                if 'islands' in pop_params.mortality.keys():
                    is_params = pop_params.mortality.islands
                    make_is = is_params.pop('make')
                    if make_is:
                        move_rast[move_rast < is_params.island_val] = 1e-8
                            #TODO: if I keep this functionality, consider adding a params option to read in the mask
                            #from a raster file or numpy array or both
                        #set the scape's mask_island_vals attribute
                        land[move_surf_scape_num].mask_island_vals = True
                            #TODO: if I keep this functionality, also come up with a better way to tell the land 
                            #to mask the 1e-8 vals when plotting the raster
                        #create a new scape for the island mask
                        is_mask = landscape.Scape(np.ma.masked_less_equal(move_rast, 1e-8).mask)
                        #set its island_mask attribute to True
                        is_mask.island_mask = True
                        #then set it as the last scape in the Land object
                        land[land.n_scapes] = is_mask
                        #and set the land's island_mask_scape_num attribute
                        land.island_mask_scape_num = land.n_scapes
                        #and increment land.n_scapes by 1
                        land.n_scapes += 1
                #reset the raster from the scape the move_surf is based on back to the raster we grabbed from
                #it earlier, in case any changes were made by the islands functionality
                    #TODO: because multiple populations could use the same scape for different things, it actually
                    #makes no sense to do this any more!
                land[move_surf_scape_num].rast = move_rast
######################################################################################################################

                #make the movement surface and set it as the pop's move_surf attribute
                pop.move_surf = spt.Movement_Surface(land[move_surf_scape_num], **ms_params)

    #if this population has changes parameterized, create a Pop_Changer object for it
    if 'change' in pop_params.keys():
        #grab the change params
        ch_params = pop_params.change
        #make Pop_Changer and set it to the pop's changer attribute
        if land.changer is not None:
            pop.changer = change.Pop_Changer(pop, ch_params, land = land)
        else:
            pop.changer = change.Pop_Changer(pop, ch_params, land = None)

    return pop 


# function for reading in a pickled pop
def read_pickled_pop(filename):
    import cPickle
    with open(filename, 'rb') as f:
        pop = cPickle.load(f)
    return pop

