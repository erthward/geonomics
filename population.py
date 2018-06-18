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

import genome
import individual
import movement
import mating
import gametogenesis
import dispersal
import selection
import mutation
import landscape
import demography

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


# ------------------------------------
# CLASSES ---------------------------
# ------------------------------------


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

        self.n_births = []  # tracker of number of births each time pop.mate is called

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

    def census(self):
        return len(self.individs)

    # method to set self.K
    def set_K(self, K):  # NOTE: Requires a landscape.Landscape instance
        self.K = K

    # method to set self.N
    def set_N(self, N):  # NOTE: Requires a landscape.Landscape instance
        self.N = N

    # method to append current pop size to the pop.Nt list
    def set_Nt(self):
        self.Nt.append(self.census())

    # method to increment all population's age by one (also adds current pop size to tracking array)
    def increment_age_stage(self, burn=False):

        # increment age of all individuals
        [ind.increment_age_stage() for ind in self.individs.values()];

        # add 1 to pop.t
        if burn == False:  # as long as this is not during burn-in, pop.t will increment
            self.t += 1

    # method to move all individuals simultaneously, and sample their new habitats
    def move(self, land):
        movement.move(self, land)
        self.set_habitat(land)
        self.set_coords_and_cells()


    # function for finding all the mating pairs in a population
    def find_mating_pairs(self, land):

        mating_pairs = mating.find_mates(self, land)
        return (mating_pairs)

    # function for executing mating for a population
    def mate(self, land, mating_pairs, burn=False):

        n_births = mating.determine_n_births(len(mating_pairs), self.n_births_lambda)
        total_births = sum(n_births)
        self.n_births.append(total_births)

        if burn == False:
            recomb_paths = self.genomic_arch.recomb_paths.get(total_births*2)
            new_genomes = mating.mate(self, mating_pairs, n_births, recomb_paths)
        
        for n_pair, pair in enumerate(mating_pairs):

            parent_centroid_x = (self.individs[pair[0]].x + self.individs[pair[1]].x)/2
            parent_centroid_y = (self.individs[pair[0]].y + self.individs[pair[1]].y)/2

            n_offspring = n_births[n_pair]
            n_gametes = 2 * n_offspring

            offspring_key = next(reversed(self.individs)) + 1
            for n in range(n_offspring):

                offspring_x, offspring_y = dispersal.disperse(land, parent_centroid_x, parent_centroid_y,
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
    def mutate(self, log=False):
        mutation.mutate(self)

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

        Valid values for normalize_by currently include 'census' and 'none'. If normalize_by = 'census', max_1 =
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

    # NOTE: DEH 01-11-18: Reformatting the habitat-getting approach
    # 1.) determined a 10x-faster way of setting hab vals of individs
    # 2.) Changed query_habitat() to set_habitat()
    # 2.) creating separate set functions for whole pop (usually called this way, so avoids conditional
    # tests) and for scape_num and/or individ-specific calls (and ditto for the get functions)


    def set_habitat(self, land, individs = None):
        if individs is None:
            inds_to_set = self.individs.values()
        
        else:
            ig = itemgetter(*individs)
            inds_to_set = ig(self.individs)
            if type(inds_to_set) is individual.Individual:
                inds_to_set = (inds_to_set,)
       
        hab = [ind.set_habitat([sc.raster[int(ind.y), int(ind.x)] for sc in list(land.scapes.values())]) for ind in inds_to_set]
        

    #method to set population's coords and cells arrays
    def set_coords_and_cells(self):
        self.coords = self.get_coords()
        self.cells = np.int8(self.coords)


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
        return selection.get_fitness(self, trait_num = trait_num, set_fit = set_fit)

    def hist_fitness(self):
        plt.hist(list(selection.get_fitness(self).values()))

    def get_single_trait_fitness(self, trait):
        return selection.get_fitness(self, trait_num = trait)

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


    def show(self, land, scape_num=None, color='black', colorbar=True, markersize=25, im_interp_method='nearest',
             alpha=False):
        # if land != None:
        if scape_num != None:
            land.scapes[scape_num].show(colorbar=colorbar, im_interp_method=im_interp_method, pop=True)
        else:
            land.show(colorbar=colorbar, im_interp_method=im_interp_method, pop=True)

        c = np.array(list(self.get_coords()))
        # NOTE: subtract 0.5 to line up the points with the plt.imshow() grid of the land; imshow plots each pixel centered on its index, but the points then plot against those indices, so wind up shifted +0.5 in each axis
        x = c[:, 0] - 0.5
        y = c[:, 1] - 0.5
        if alpha == True:
            alpha = 0.6
        else:
            alpha = 1.0

        plt.scatter(x, y, s=markersize, c=color, alpha=alpha);
        plt.xlim(-0.6, land.dims[1] - 0.4);
        plt.ylim(-0.6, land.dims[0] - 0.4)

    # mpl.pyplot.plot([n[0] for n in coords], [n[1] for n in coords], 'ko', scalex = False, scaley = False, color = color, markersize = markersize, alpha = alpha)

    def show_individs(self, individs, land, scape_num=None, color='black', im_interp_method='nearest', markersize=40,
                      alpha=0.5):
        # if land != None and scape_num != None:
        land.scapes[scape_num].show(im_interp_method=im_interp_method, pop=True)

        # coords = dict([(k, (ind.x, ind.y)) for k, ind in self.individs.items() if k in individs])
        c = np.array(list(self.get_coords(individs)))
        # NOTE: subtract 0.5 to line up points with imshow grid; see note in the pop.show() definition for details
        x = c[:, 0] - 0.5
        y = c[:, 1] - 0.5
        plt.scatter(x, y, s=markersize, c=color, alpha=alpha);
        plt.xlim(-0.6, land.dims[1] - 0.4);
        plt.ylim(-0.6, land.dims[0] - 0.4)

    # NOTE: perhaps worth figuring out how to label with the individual number!!
    # for k, coord_pair in coords.items():
    # ax = mpl.pyplot.plot(coord_pair[0], coord_pair[1], 'ko', scalex = False, scaley = False, color = color, markersize = 8.5)

    def show_density(self, land, normalize_by='census', max_1=False, color='black', markersize=40, alpha=0.5):
        dens = self.calc_density(land, normalize_by=normalize_by, max_1=max_1)
        dens.show(im_interp_method='nearest', pop=True)

        c = np.array(list(self.get_coords()))
        # NOTE: subtract 0.5 to line up points with imshow grid; see note in the pop.show() definition for details
        x = c[:, 0] - 0.5
        y = c[:, 1] - 0.5
        plt.scatter(x, y, s=markersize, c=color, alpha=alpha);
        plt.xlim(-0.6, land.dims[1] - 0.4);
        plt.ylim(-0.6, land.dims[0] - 0.4)

    # ax = mpl.pyplot.plot([i[0] - 0.5 for i in list(c.values())], [i[1] - 0.5 for i in list(c.values())], 'ko', scalex = False, scaley = False, color = color, markersize = 8.5)

    # NOTE: perhaps worth figuring out how to label with the individual number!!

    # method for plotting individuals colored by their genotype at a given locus
    def show_genotype(self, locus, land, scape_num=None, im_interp_method='nearest', markersize=65, alpha=1,
                      by_dominance=False):
        if scape_num != None:
            land.scapes[scape_num].show(im_interp_method=im_interp_method, pop=True)

        else:
            land.show(im_interp_method=im_interp_method, pop=True)

        if by_dominance == True:
            genotypes = self.get_genotype(locus, by_dominance=True)
        else:
            genotypes = self.get_genotype(locus)

        colors = ['#3C22B4', '#80A6FF',
                  '#FFFFFF']  # COLORS TO MATCH LANDSCAPE PALETTE EXTREMES, BUT WITH HYBRID A MIX OF THE EXTREMES RATHER THAN THE YELLOW AT THE MIDDLE OF THE PALETTE, FOR NICER VIEWING: blue = [0,0], light blue = [0,1], white = [1,1]
        # colors = ['#ff4d4d', '#ac72ac', '#4d4dff'] # red = [0,0], purple = [0,1], blue = [1,1]
        for n, genotype in enumerate([0.0, 0.5, 1.0]):
            inds = [i for i, g in genotypes.items() if np.atleast_1d(g)[0] == genotype]
            # plot if there are any individuals of this genotype
            if len(inds) >= 1:
                c = self.get_coords(inds)
                x = c[:, 0] - 0.5
                y = c[:, 1] - 0.5
                plt.scatter(x, y, s=markersize, c=colors[n], alpha=alpha);
                plt.xlim(-0.6, land.dims[1] - 0.4);
                plt.ylim(-0.6, land.dims[0] - 0.4)

    # method for plotting individuals colored by their phenotypes for a given trait
    def show_phenotype(self, trait, land, scape_num=None, im_interp_method='nearest', markersize=65, alpha=1):

        if scape_num != None:
            land.scapes[scape_num].show(im_interp_method=im_interp_method, pop=True)
        else:
            land.scapes[self.genomic_arch.traits[trait].scape_num].show(im_interp_method=im_interp_method, pop=True)

        from matplotlib.colors import LinearSegmentedColormap
        colors = ['#3C22B4', '#80A6FF', '#FFFFFF']
        # colors to match landscape palette extremes, but with hybrid a mix of the extremes
        # rather than the yellow at the middle of the palette, for nicer viewing: blue = [0,0],
        # light blue = [0,1], white = [1,1]
        cmap = LinearSegmentedColormap.from_list('my_cmap', colors, N=50)

        c = self.get_coords()
        z = {i:v[trait] for i,v in self.get_phenotype().items()}

        data = list(OD({i: (c[i][0] - 0.5, c[i][1] - 0.5, z[i]) for i in self.individs.keys()}).values())

        plt.scatter([i[0] for i in data], [i[1] for i in data], s=markersize, c=[i[2] for i in data], cmap=cmap,
                    linewidth=1, edgecolor='black',
                    alpha=alpha)
        plt.xlim(-0.6, land.dims[1] - 0.4)
        plt.ylim(-0.6, land.dims[0] - 0.4)
        # plt.scatter(x,y, s = markersize, c = [z[i] for i in i.items()}nds], cmap = cmap, alpha = alpha);plt.xlim(-0.6, land.dims[1]-0.4); plt.ylim(-0.6, land.dims[0]-0.4)

    # method for plotting individuals colored and sized by their overall fitnesses
    def show_fitness(self, land, scape_num=None, im_interp_method='nearest', min_markersize=60, alpha=1):

        if scape_num != None:
            land.scapes[scape_num].show(im_interp_method=im_interp_method, pop=True)
        else:
            land.scapes[land.n_movement_surf_scape].show(im_interp_method=im_interp_method, pop=True)

        from matplotlib.colors import LinearSegmentedColormap
        colors = ['#3C22B4', '#80A6FF', '#FFFFFF']
        # colors to match landscape palette extremes, but with hybrid a mix of the extremes
        # rather than the yellow at the middle of the palette, for nicer viewing: blue = [0,0],
        # light blue = [0,1], white = [1,1]
        cmap = LinearSegmentedColormap.from_list('my_cmap', colors, N=50)

        # get all individs' coords
        c = self.get_coords()

        # get all individs' fitness values
        w = self.get_fitness()

        # calc minimum possible fitness
        min_fit = np.product([1 - np.atleast_2d(t.phi).min() for t in list(self.genomic_arch.traits.values())])

        # generate an evenly spaced range of the possible fitness values
        fit_vals = np.linspace(min_fit, 1, 50)

        # use index of closest possible fitness val to get a markersize differential (to be added to min markersize) for each individual
        markersize_differentials = {i: 3 * np.abs(fit_vals - w[n]).argmin() for n,i in enumerate(self.individs.keys())}

        data = list(OD({i: (c[i][0] - 0.5, c[i][1] - 0.5, w[i], markersize_differentials[i]) for i in
                        self.individs.keys()}).values())

        plt.scatter([i[0] for i in data], [i[1] for i in data], s=[min_markersize + i[3] for i in data],
                    c=[i[2] for i in data], cmap=cmap, alpha=alpha)
        plt.xlim(-0.6, land.dims[1] - 0.4)
        plt.ylim(-0.6, land.dims[0] - 0.4)

        # plt.scatter(x,y, s = [min_markersize + diff for diff in markersize_differentials] , c = [fits[i] for i in inds], cmap = cmap, alpha = alpha);plt.xlim(-0.6, land.dims[1]-0.4); plt.ylim(-0.6, land.dims[0]-0.4)

    # method for plotting individuals colored by their phenotypes for a given trait, sized by their fitness
    def show_single_trait_fitness(self, trait, land, scape_num=None, im_interp_method='nearest', min_markersize=60,
                                  alpha=1):

        if scape_num != None:
            land.scapes[scape_num].show(im_interp_method=im_interp_method, pop=True)
        else:
            land.scapes[self.genomic_arch.traits[trait].scape_num].show(im_interp_method=im_interp_method, pop=True)

        from matplotlib.colors import LinearSegmentedColormap
        colors = ['#3C22B4', '#80A6FF', '#FFFFFF']
        # colors to match landscape palette extremes, but with hybrid a mix of the extremes
        # rather than the yellow at the middle of the palette, for nicer viewing: blue = [0,0],
        # light blue = [0,1], white = [1,1]
        cmap = LinearSegmentedColormap.from_list('my_cmap', colors, N=50)

        # get all individs' coords
        c = self.get_coords()

        # calc minimum possible fitness
        min_fit = 1 - np.atleast_2d(self.genomic_arch.traits[trait].phi).min()

        # generate an evenly spaced range of the possible fitness values
        fit_vals = np.linspace(min_fit, 1, 50)

        # get all individs' fitness values
        w = self.get_single_trait_fitness(trait)

        # use index of closest possible fitness val to get a markersize differential (to be added to min markersize) for each individual
        markersize_differentials = [3 * np.abs(fit_vals - w[i]).argmin() for i in range(self.census())]

        z = [v[trait] for v in self.get_phenotype()]

        data = list(OD({i: (c[n][0] - 0.5, c[n][1] - 0.5, w[n], z[n], markersize_differentials[n]) for n,i in
                        enumerate(self.individs.keys())}).values())

        # separate colormap to color marker edges from black (fit = 1) to white (fit = 0) through red
        inside_marker_cmap = mpl.cm.get_cmap('RdYlGn')

        plt.scatter([i[0] for i in data], [i[1] for i in data], s=[min_markersize + i[4] for i in data],
                    c=[i[3] for i in data], cmap=cmap, alpha=alpha)
        plt.scatter([i[0] for i in data], [i[1] for i in data], s=min_markersize / 3, c=[i[2] for i in data],
                    cmap=inside_marker_cmap)
        plt.xlim(-0.6, land.dims[1] - 0.4)
        plt.ylim(-0.6, land.dims[0] - 0.4)

        # plt.scatter(x,y, s = [min_markersize + diff for diff in markersize_differentials] , c = [z[i] for i in inds], cmap = cmap, alpha = alpha)
        # plt.scatter(x,y, s = min_markersize/3 , c = [fits[i] for i in inds], cmap = inside_marker_cmap)

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

    # method to plot (or add to an open plot) individuals' IDs
    def show_ind_ids(self):
        [plt.text(v[0] - 0.5, v[1] - 0.5, i) for i, v in dict(zip(self.individs.keys(), self.get_coords()))]

    def pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)


# --------------------------------------
# FUNCTIONS ---------------------------
# --------------------------------------


def create_population(genomic_arch, land, params, burn=False):
    # grab necessary params from params dict

    dims = land.dims

    assert dims.__class__.__name__ in ['tuple', 'list'], "dims should be expressed as a tuple or a list"
    individs = OD()
    
    N = params['pop']['main']['N_start']

    for idx in range(N):
        # use individual.create_individual to simulate individuals and add them to the population
        ind = individual.create_individual(idx = idx, genomic_arch = genomic_arch, dims = dims)
        individs[idx] = ind
        #land.mating_grid.add(ind)

    pop = Population(p_params = params['pop'],  individs=individs, genomic_arch=genomic_arch)

    #set initial habitat values
    pop.set_habitat(land)

    #set initial coords and cells
    pop.set_coords_and_cells()

    #set phenotypes
    [i.set_phenotype(pop.genomic_arch) for i in pop.individs.values()]

    return pop


# function for reading in a pickled pop
def load_pickled_pop(filename):
    import cPickle
    with open(filename, 'rb') as f:
        pop = cPickle.load(f)

    return pop
