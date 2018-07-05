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
import viz
import spatial as spt
import genome
import individual
import movement
import mating
import selection
import mutation
import landscape
import demography

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

class Population:
    '''Population(self, N, individs, genomic_arch, size, birth_rate = None, death_rate = None, T = None)'''

    def __init__(self, p_params, individs, genomic_arch):

        self.it = None  # attribute to keep track of iteration number this pop is being used for
        # (optional; will be set by the iteration.py module if called)

        self.t = 0  # attribute to keep of track of number of timesteps the population has evolved
        # NOTE: This way, if model is stopped, changes are made, then it is run further,
        # this makes it easy to continue tracking from the beginning

        #grab all of the params['pop'] parameters as Population attributes
        for section in ['main', 'mating', 'mortality', 'movement']:
            for att in p_params[section].keys():
                setattr(self, att, p_params[section][att])

        self.N = None  # attribute to hold a landscape.Landscape object of the current population density

        self.K = None  # attribute to hold an landscape.Landscape object of the local carrying capacity (i.e. 'target' dynamic equilibrium population density)

        self.Nt = []  # list to record population size (appended each time
        # pop.increment_age_stage() is called)

        self.n_births = []  # tracker of number of births each time pop.do_mating is called

        self.n_deaths = []  # tracker of number of deaths each time demography.pop_dynamics is called

        self.individs = individs  # dict of all individuals, as instances of individual.Individual class

        self.initial_size = len(self.individs)

        #attributes for storing numpy arrays of all individuals' coordinates and cells, to avoid repeat compute time each turn
        self.coords = None
        self.cells = None

        self.genomic_arch = genomic_arch
        # Add other demographic parameters?? Other stuff in general??

        self.heterozygote_effects = {
            0: lambda g: [np.ceil(np.mean(g[i,])) for i in range(g.shape[0])],
            # relative fitness of homozygous 1 = relative fitness of heterozygote, i.e. 1 = 1-hs
            0.5: lambda g: [np.mean(g[i,]) for i in range(g.shape[0])],
            # relative fitness of heterozygote halfway between homozygotes, i.e 1-hs = 1-s/2
            1: lambda g: [np.floor(np.mean(g[i,])) for i in range(g.shape[0])]
            # relative fitness of homozygous 0 = relative fitness of heterozygote, i.e. 1-hs = 1-s
        }

        #create empty attribute to hold kd_tree
        self.kd_tree = None

        assert type(self.N_start) == int, "N must be an integer"
        assert type(self.individs) == OD, "self.individs must be of type collections.OrderedDict"
        assert list(set([i.__class__.__name__ for i in list(self.individs.values())])) == [
            'Individual'], "self.individs must be a dictionary containing only instances of the individual.Individual class"
        assert self.genomic_arch.__class__.__name__ == 'Genomic_Architecture', "self.genomic_arch must be an instance of the genome.Genomic_Architecture class"


        self.__coord_attrgetter__ = attrgetter('x', 'y')
        self.__individ_coord_attrgetter__ = attrgetter('idx', 'x', 'y')

    #####################
    ### OTHER METHODS ###
    #####################

    def get_size(self):
        return len(self.individs)

    # method to set self.K
    def set_K(self, K):  # NOTE: Requires a landscape.Landscape instance
        self.K = K

    # method to set self.N
    def set_N(self, N):  # NOTE: Requires a landscape.Landscape instance
        self.N = N

    # method to append current pop size to the pop.Nt list
    def set_Nt(self):
        self.Nt.append(self.get_size())

    # method to increment all population's age by one (also adds current pop size to tracking array)
    def reset_age_stage(self, burn=False):

        # increment age of all individuals
        [ind.reset_age_stage() for ind in self.individs.values()];

        # add 1 to pop.t
        if burn == False:  # as long as this is not during burn-in, pop.t will increment
            self.t += 1

    # method to move all individuals simultaneously, and sample their new habitats
    def do_movement(self, land):
        movement.move(self, land)
        self.set_habitat(land)
        self.set_coords_and_cells()


    # function for finding all the mating pairs in a population
    def find_mating_pairs(self, land):

        mating_pairs = mating.find_mates(self, land)
        return (mating_pairs)

    # function for executing mating for a population
    def do_mating(self, land, mating_pairs, burn=False):

        n_births = mating.draw_n_births(len(mating_pairs), self.n_births_lambda)
        total_births = sum(n_births)
        self.n_births.append(total_births)

        if burn == False:
            recomb_paths = self.genomic_arch.recomb_paths.get_paths(total_births*2)
            new_genomes = mating.do_mating(self, mating_pairs, n_births, recomb_paths)
        
        for n_pair, pair in enumerate(mating_pairs):

            parent_centroid_x = (self.individs[pair[0]].x + self.individs[pair[1]].x)/2
            parent_centroid_y = (self.individs[pair[0]].y + self.individs[pair[1]].y)/2

            n_offspring = n_births[n_pair]
            n_gametes = 2 * n_offspring

            offspring_key = next(reversed(self.individs)) + 1
            for n in range(n_offspring):

                offspring_x, offspring_y = movement.disperse(land, parent_centroid_x, parent_centroid_y,
                        self.dispersal_mu, self.dispersal_sigma)

                if self.sex:
                    offspring_sex = r.binomial(1, 0.5)

                age = 0

                if burn == False:
                    if self.sex == True:
                        self.individs[offspring_key] = individual.Individual(offspring_key, new_genomes[n_pair][n], offspring_x, offspring_y, offspring_sex, age)
                    else:
                        self.individs[offspring_key] = individual.Individual(offspring_key, new_genomes[n_pair][n], offspring_x, offspring_y, age)
               
                    #set new individual's phenotype (won't be set during burn-in, becuase no genomes assigned)
                    self.individs[offspring_key].set_phenotype(self.genomic_arch)


                elif burn:
                    if self.sex == True:
                        self.individs[offspring_key] = individual.Individual(offspring_key, np.array([0]), offspring_x, offspring_y, offspring_sex, age)
                    else:
                        self.individs[offspring_key] = individual.Individual(offspring_key, np.array([0]), offspring_x, offspring_y,
                                                                             age)
                #land.mating_grid.add(self.individs[offspring_key])

                offspring_key += 1

                
        # sample all individuals' habitat values, to initiate for offspring
        self.set_habitat(land)
        self.set_coords_and_cells()

        print('\n\t%i individuals born' % (total_births))

    # method to carry out mutation
    def do_mutation(self, log=False):
        mutation.do_mutation(self)

    def check_extinct(self):
        if len(self.individs.keys()) == 0:
            print('\n\nYOUR POPULATION WENT EXTINCT!\n\n')
            return (1)
            # sys.exit()
        else:
            return (0)

    def calc_density(self, land, normalize_by= None, min_0=True, max_1=False, max_val=None, as_landscape = False, set_N=False):

        '''
        Calculate an interpolated raster of local population density, using the
        Landscape_Stack.density_grid_stack object.

        Valid values for normalize_by currently include 'pop_size' and 'none'. If normalize_by = 'pop_size', max_1 =
        True will cause the output density raster to vary between 0 and 1, rather than between 0 and the current
        max normalized density value. 
        '''

        x = self.get_x_coords()
        y = self.get_y_coords()

        dens = land.density_grid_stack.calc_density(x, y)

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
            dens = landscape.Landscape(land.dims, dens)

        if set_N:
            self.set_N(dens)

        else:
            return (dens)


    def set_habitat(self, land, individs = None):
        if individs is None:
            inds_to_set = self.individs.values()
        
        else:
            ig = itemgetter(*individs)
            inds_to_set = ig(self.individs)
            if type(inds_to_set) is individual.Individual:
                inds_to_set = (inds_to_set,)
       
        hab = [ind.set_habitat([sc.raster[int(ind.y), int(ind.x)] for sc in list(land.scapes.values())]) for ind in inds_to_set]

    
    def set_phenotype(self):
        [ind.set_phenotype(self.genomic_arch) for ind in self.individs.values()];


    def set_fitness(self):
        indices = list(self.individs.keys())
        fit = self.get_fitness()
        [self.individs[ind].set_fitness(fit[n]) for n, ind in enumerate(indices)]
        

    #method to set population's coords and cells arrays
    def set_coords_and_cells(self):
        self.coords = self.get_coords()
        self.cells = np.int8(self.coords)

    #method to set the population's kd_tree attribute
    def set_kd_tree(self, leafsize = 100):
        self.kd_tree = spt.KD_Tree(coords = self.coords, leafsize = leafsize)


    # method to get individs' habitat values
    def get_habitat(self, scape_num=None, individs=None):
        if individs is None:
            if scape_num is None:
                habs = np.array([ind.habitat for i, ind in self.individs.items()])
            else:
                habs = np.array([ind.habitat[scape_num] for i, ind in self.individs.items()])
        else:
            ig = itemgetter(*individs)
            if scape_num is None:
                habs = {i:ind.habitat for i, ind in self.individs.items()}
                habs = np.array(ig(habs))
            else:
                habs = {i:ind.habitat[scape_num] for i, ind in self.individs.items()}
                habs = np.array(ig(habs))
        return(habs)

    def get_age(self, individs=None):
        if individs is None:
            ages = np.array([ind.age for ind in self.individs.values()])
        else:
            ig = itemgetter(*individs)
            ages = {i: ind.age for i, ind in self.individs.items()}
            ages = np.array(ig(ages))
        return(ages)

    def get_genotype(self, locus, return_format='mean', individs=None, by_dominance=False):

        if individs is None:
            individs = self.individs.keys()
            # individs = range(len(self.genomic_arch.s[chromosome]))

        if return_format == 'biallelic':
            return {i: self.individs[i].genome[locus, :] for i in self.individs.keys() if i in individs}

        elif return_format == 'mean':
            if by_dominance == True:
                h = self.genomic_arch.h[locus]
            else:
                h = 0.5
            return dict(zip(individs, self.heterozygote_effects[h](
                np.array([ind.genome[locus,] for i, ind in self.individs.items() if i in individs]))))

    # convenience method for getting whole population's phenotype
    def get_phenotype(self, individs=None):
        if individs is None:
            zs = np.array([ind.phenotype for ind in self.individs.values()])
        else: 
            ig = itemgetter(*individs)
            zs = {i:ind.phenotype for i, ind in self.individs.items()}
            zs = np.array(ig(zs))
        return(zs)

    def get_fitness(self, trait_num = None, set_fit = False):
        return selection.calc_fitness(self, trait_num = trait_num, set_fit = set_fit)

    def show_hist_fitness(self):
        plt.hist(list(selection.get_fitness(self).values()))

    def get_dom(self, locus):
        return {locus: self.genomic_arch.h[locus]}

    def get_coords(self, individs = None, float = True):
        coords = list(map(self.__coord_attrgetter__, self.individs.values()))
        if individs is not None: 
            ig = itemgetter(*individs)
            coords = ig(dict(zip(self.individs.keys(), coords)))

        if float == True:
            coords = np.float32(coords)
        else:
            coords = np.int8(coords)

        return(coords)
            
            
    def get_cells(self, individs = None):
        cells = self.get_coords(individs = individs, float = False)
        return(cells)
  

    def get_x_coords(self, individs=None):
        coords = self.get_coords(individs = individs)
        return(coords[:,0])


    def get_y_coords(self, individs=None):
        coords = self.get_coords(individs = individs)
        return(coords[:,1])

    #use the kd_tree attribute to find nearest neighbors either within the population, if within == True,
    #or between the population and the points provided, if within == False and points is not None
    def find_neighbors(self, dist, within=True, coords=None, k = 2):
        if within:
            coords = self.coords
        dists, pairs = self.kd_tree.query(coords = coords, dist = dist, k = k)
        return(dists, pairs)

    #method for plotting the population (or a subset of its individuals, by ID) on top of a landscape (or landscape stack)
    def show(self, land, scape_num=None, hide_land=False, individs=None, text=False, color='black',
            edge_color='face', text_color='black', colorbar=True, size=25, text_size=9, im_interp_method='nearest', 
            land_cmap='terrain', pt_cmap=None, alpha=False, zoom_width=None, x=None, y=None):
        #get coords
        if individs is None:
            coords = self.coords
            if text:
                text = list(self.individs.keys())
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
        viz.show_points(coords, scape_num = scape_num, color = color, edge_color = edge_color, text_color = text_color, size = size, text_size = text_size, 
                alpha = alpha, text = text, plt_lims = plt_lims, pt_cmap = pt_cmap)


    #method for plotting the population on top of its estimated population-density raster
    def show_density(self, land, normalize_by='pop_size', individs=None, text=False, max_1=False, color='black', edge_color='face', 
            text_color='black', size=25, text_size = 9, alpha=0.5, zoom_width=None, x=None, y=None):
        dens = self.calc_density(land, normalize_by=normalize_by, max_1=max_1)
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
                    im_interp_method = im_interp_method, alpha = alpha, zoom_width = zoom_width, x = x, y = y)


    # method for plotting individuals colored by their phenotypes for a given trait
    def show_phenotype(self, trait, land, scape_num=None, individs=None, text=False, size=25, text_size = 9, 
            edge_color='black', text_color='black', colorbar=True, im_interp_method='nearest', 
            alpha=1, by_dominance=False, zoom_width=None, x=None, y=None):

        z = dict(zip(self.individs.keys(), self.get_phenotype()[:,trait]))
        if individs is not None:
            z = {i:v for i,v in z.items() if i in individs}

        self.show(land, scape_num = scape_num, individs = individs, text = text, color = list(z.values()),
                pt_cmap = 'terrain', edge_color = edge_color, text_color = text_color, colorbar = colorbar, 
                size = size, text_size = text_size, im_interp_method = im_interp_method, alpha = alpha, zoom_width = zoom_width, x = x, y = y)


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
        w = dict(zip(self.individs.keys(), w))
        if individs is not None:
            w = {i:v for i,v in w.items() if i in individs}

        # calc minimum possible fitness (for phenotypes within 0 <= z <= 1, which in reality isn't a
        # constraint, but values lower than this will also be constrained to the minimum-value color for plotting)
        if trait_num is None: 
            min_fit = np.product([1 - np.atleast_2d(t.phi).min() for t in list(self.genomic_arch.traits.values())])
        else:
            min_fit = 1 - np.atleast_2d(self.genomic_arch.traits[trait_num].phi).min()

        #then get uneven cmap and cbar-maker (needs to be uneven to give color-resolution to values varying
        #between 1 and the minimum-fitness value assuming all phenotypes are constrained 0<=z<=1, but then also
        #allow gradually deepening reds for fitness values lower than that minimum value), using the min_fit val
        cmap, make_cbar = viz.make_fitness_cmap_and_cbar_maker(min_val = min_fit, max_val = 1, cmap = fit_cmap, trait_num = trait_num)
        
        #plot the trait phenotype in larger circles first, if trait is not None
        if trait_num is not None:
            self.show_phenotype(trait=trait_num, land=land, scape_num=scape_num, individs=individs, text=False,
                    size=size, text_size = text_size, edge_color=edge_color, text_color=text_color,
                    colorbar=colorbar, im_interp_method=im_interp_method, 
                    alpha=alpha, by_dominance=by_dominance, zoom_width=zoom_width, x=x, y=y)
            size = round(0.2*size)

        self.show(land, scape_num = scape_num, individs = individs, text = text, color = list(w.values()),
                pt_cmap = cmap, edge_color = edge_color, text_color = text_color, colorbar = colorbar, 
                size = size, text_size = text_size, im_interp_method = im_interp_method, alpha = alpha, zoom_width = zoom_width, x = x, y = y)

        #and make a colorbar for the fitness values 
        viz.make_fitness_cbar(make_cbar, min_fit)


    # method for plotting a population pyramid
    # NOTE: NEED TO FIX THIS SO THAT EACH HALF PLOTS ON OPPOSITE SIDES OF THE Y-AXIS
    def show_pyramid(self):
        plt.hist([ind.age for ind in list(self.individs.values()) if ind.sex == 0], orientation='horizontal',
                 color='pink',
                 alpha=0.6)
        plt.hist([ind.age for ind in list(self.individs.values()) if ind.sex == 1], orientation='horizontal',
                 color='skyblue',
                 alpha=0.6)


    def show_pop_growth(self, params):
        T = range(len(self.Nt))
        x0 = self.N_start / self.K.sum()
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

def make_population(genomic_arch, land, params, burn=False):
    # grab necessary params from params dict

    dims = land.dims

    assert dims.__class__.__name__ in ['tuple', 'list'], "dims should be expressed as a tuple or a list"
    individs = OD()
    
    N = params['pop']['main']['N_start']

    for idx in range(N):
        # use individual.create_individual to simulate individuals and add them to the population
        ind = individual.make_individual(idx = idx, genomic_arch = genomic_arch, dims = dims)
        individs[idx] = ind
        #land.mating_grid.add(ind)

    pop = Population(p_params = params['pop'],  individs=individs, genomic_arch=genomic_arch)

    #set initial habitat values
    pop.set_habitat(land)

    #set initial coords and cells
    pop.set_coords_and_cells()

    #set phenotypes
    [i.set_phenotype(pop.genomic_arch) for i in pop.individs.values()]

    #set the kd_tree
    pop.set_kd_tree()

    return pop


# function for reading in a pickled pop
def read_pickled_pop(filename):
    import cPickle
    with open(filename, 'rb') as f:
        pop = cPickle.load(f)

    return pop
