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
import mating
import gametogenesis
import dispersal
import selection
import mutation
import landscape
import demography

import numpy as np
from numpy import random as r
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import cKDTree
from collections import Counter as C

import sys


# ------------------------------------
# CLASSES ---------------------------
# ------------------------------------


class Population:
    '''Population(self, N, individs, genomic_arch, size, birth_rate = None, death_rate = None, T = None)'''

    def __init__(self, N, individs, genomic_arch, size, T):

        self.T = T

        self.it = None  # attribute to keep track of iteration number this pop is being used for
        # (optional; will be set by the iteration.py module if called)

        self.t = 0  # attribute to keep of track of number of timesteps the population has evolved
        # NOTE: This way, if model is stopped, changes are made, then it is run further,
        # this makes it easy to continue tracking from the beginning

        self.N = None  # attribute to hold a landscape.Landscape object of the current population density

        self.K = None  # attribute to hold an landscape.Landscape object of the local carrying capacity (i.e. 'target' dynamic equilibrium population density)

        self.Nt = []  # list to record population size (appended each time
        # pop.increment_age_stage() is called)

        self.n_births = []  # tracker of number of births each time pop.mate is called

        self.n_deaths = []  # tracker of number of deaths each time demography.pop_dynamics is called

        self.individs = individs  # dict of all individuals, as instances of individual.Individual class

        self.initial_size = len(self.individs)

        if size.__class__.__name__ in ['float', 'int'] and T != None:
            self.size = [size] * T
        elif size.__class__.__name__ in ['list']:
            assert T != None and len(
                size) == T, "if expressing population size as a list, total model runtime must be provided and must equal the list length"
            self.size = size
        else:
            self.size = size

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

        assert type(N) == int, "N must be an integer"
        assert type(self.individs) == dict, "self.individs must be a dictionary"
        assert list(set([i.__class__.__name__ for i in self.individs.values()])) == [
            'Individual'], "self.individs must be a dictionary containing only instances of the individual.Individual class"
        assert self.genomic_arch.__class__.__name__ == 'Genomic_Architecture', "self.genomic_arch must be an instance of the genome.Genomic_Architecture class"

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

    # method to move all individuals simultaneously
    def move(self, land, params):
        [ind.move(land, params) for ind in self.individs.values()];

        self.set_habitat(land)

    # function for finding all the mating pairs in a population
    def find_mating_pairs(self, land, params):

        mating_pairs = mating.find_mates(self, params)
        return (mating_pairs)

    # function for executing mating for a population
    def mate(self, land, params, mating_pairs, burn=False):

        # pull necessary parameters from params dict
        mu_dispersal = params['mu_dispersal']
        sigma_dispersal = params['sigma_dispersal']
        sex = params['sex']
        repro_age = params['repro_age']

        num_births = mating.determine_num_births(len(mating_pairs), params)
        total_births = sum(num_births)
        self.n_births.append(total_births)

        if burn == False:
            recombinants = gametogenesis.recombine(self.genomic_arch.r_lookup, 2 * total_births)

        for pair in mating_pairs:

            parent_centroid_x = np.mean((self.individs[pair[0]].x, self.individs[pair[1]].x))
            parent_centroid_y = np.mean((self.individs[pair[0]].y, self.individs[pair[1]].y))

            n_offspring = num_births.pop()
            n_gametes = 2 * n_offspring

            # gamete_recomb_paths, recombinants = recombinants[:,0:n_gametes], recombinants[:,n_gametes:]
            if burn == False:
                gamete_recomb_paths, recombinants = [i.flatten() for i in
                                                     np.hsplit(recombinants[:, 0:n_gametes], n_gametes)], recombinants[
                                                                                                          :, n_gametes:]
                new_genomes = mating.mate(self, pair, params, n_offspring, gamete_recomb_paths)

            for n in range(n_offspring):

                offspring_x, offspring_y = dispersal.disperse(land, parent_centroid_x, parent_centroid_y, mu_dispersal,
                                                              sigma_dispersal)

                if sex == True:
                    offspring_sex = r.binomial(1, 0.5)

                age = 0

                offspring_key = max(self.individs.keys()) + 1

                if burn == False:
                    if sex == True:
                        ind = individual.Individual(new_genomes[n], offspring_x, offspring_y,
                                                                             offspring_sex, age)
                    else:
                        ind = individual.Individual(new_genomes[n], offspring_x, offspring_y,
                                                                             age)

                elif burn == True:
                    if sex == True:
                        ind = individual.Individual(np.array([0]), offspring_x, offspring_y,
                                                                             offspring_sex, age)
                    else:
                        ind = individual.Individual(np.array([0]), offspring_x, offspring_y,
                                                                             age)
                self.individs[offspring_key] = ind
                land.mating_grid.add(ind)

        # sample all individuals' habitat values, to initiate for offspring
        self.set_habitat(land)

        print '\n\t%i individuals born' % (total_births)

    # method to carry out mutation
    def mutate(self, log=False):
        mutation.mutate(self)

    def check_extinct(self):
        if len(self.individs.keys()) == 0:
            print '\n\nYOUR POPULATION WENT EXTINCT!\n\n'
            return (1)
            # sys.exit()
        else:
            return (0)

    def calc_density(self, land, window_width=None, normalize_by='none', min_0=True, max_1=False, max_val=None,
                     set_N=False):

        '''
        Calculate an interpolated raster of local population density, using a window size of window_width.
        Valid values for normalize_by currently include 'census' and 'none'. If normalize_by = 'census', max_1 =
        True will cause the output density raster to vary between 0 and 1, rather than between 0 and the current
        max normalized density value. Window width will default to 1/10 of the larger raster dimension.
        '''

        # window width defaults to 1/10 the maximum landscape dimension
        if window_width == None:
            window_width = max(land.dims) * 0.1

        # shorthand
        dims = land.dims

        # get a list of pop's coord-tuples
        c = self.get_coords().values()

        # make window_width a float, to avoid Py2 integer-division issues
        window_width = float(window_width)

        # create meshgrid using window_width/2 as step size
        grid_j, grid_i = np.mgrid[0:dims[0]:complex("%ij" % (dims[0] / (window_width / 2))),
                         0:dims[1]:complex("%ij" % (dims[1] / (window_width / 2)))]

        # grid_j, grid_i = np.mgrid[0+(window_width/2):dims[0]-(window_width/2):complex("%ij" % (dims[0]/(window_width/2))), 0+(window_width/2):dims[1]-(window_width/2):complex("%ij" % (dims[1]/(window_width/2)))]

        # flatten the arrays, so that I can run over them in a single for loop
        gj = grid_j.ravel()
        gi = grid_i.ravel()

        # make lists of tuples, of same length as gj, containing the window ll and ur coords
        window_ll = [(max(gj[n] - (window_width / 2), 0), max(gi[n] - (window_width / 2), 0)) for n in
                     range(len(gj))]  # constrain min window vals to 0
        window_ur = [(min(gj[n] + (window_width / 2), land.dims[0]), min(gi[n] + (window_width / 2), land.dims[1])) for
                     n in range(len(gj))]  # constrain max window vals to each respective land dimension
        assert len(window_ll) == len(gj)
        assert len(window_ur) == len(gj)

        # make a list of the sizes of each window
        window_size = [(window_ur[n][0] - window_ll[n][0]) * (window_ur[n][1] - window_ll[n][1]) for n in
                       range(len(gj))]  # calculate size of this window (not always the same because of edge effects
        assert len(window_size) == len(gj)

        # make a list of the counts of all individs within each window
        window_ct = [len([ind for ind in range(len(c)) if
                          (c[ind][0] > window_ll[n][0] and c[ind][0] <= window_ur[n][0]) and (
                                      c[ind][1] > window_ll[n][1] and c[ind][1] <= window_ur[n][1])]) for n in
                     range(len(gj))]
        assert len(window_ct) == len(gj)

        # divide window counts by window sizes
        window_dens = [window_ct[n] / window_size[n] for n in range(len(window_ct))]  # divide by window size
        assert len(window_dens) == len(gj)

        # if normalize_by == census, then divide each density by total pop census size
        if normalize_by == 'census':
            N = self.census()
            window_dens = [dens / N for dens in window_dens]
        elif normalize_by == 'none':
            pass

        else:  # POTENTIALLY ADD OTHER OPTIONS HERE, i.e. to normalize by starting size?
            pass

            # interpolate resulting density vals to a grid equal in size to the landscape
        new_gj, new_gi = np.mgrid[0:dims[0] - 1:complex("%ij" % (dims[0])), 0:dims[1] - 1:complex("%ij" % (dims[1]))]
        dens = interpolate.griddata(np.array(zip(list(gi), list(gj))), window_dens, (new_gj, new_gi), method='cubic')

        if normalize_by != 'none':

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

        if min_0 == True:
            dens[dens < 0] = 0

        if max_val != None:
            dens[dens > max_val] = max_val

        if set_N == True:
            self.set_N(landscape.Landscape(dims, dens))

        else:
            return (landscape.Landscape(dims, dens))

    # NOTE: DEH 01-11-18: Reformatting the habitat-getting approach
    # 1.) determined a 10x-faster way of setting hab vals of individs
    # 2.) Changed query_habitat() to set_habitat()
    # 2.) creating separate set functions for whole pop (usually called this way, so avoids conditional
    # tests) and for scape_num and/or individ-specific calls (and ditto for the get functions)

    # method for setting all individs' habitat attributes to match current locations
    def set_habitat(self, land, scape_num=None, individs=None):
        hab = {i: [sc.raster[int(ind.y), int(ind.x)] for sc in land.scapes.values()] for i, ind in
               self.individs.items()}
        for i, v in self.individs.items():
            v.habitat = hab[i]

    # if scape_num and individs args set to None, does same as set_habitat()
    def set_habitat_by_land_ind(self, land, scape_num=None, individs=None):
        if individs is None:
            if scape_num is None:
                hab = {i: [sc.raster[int(ind.y), int(ind.x)] for sc in land.scapes.values()] for i, ind in
                       self.individs.items()}
            else:
                hab = {i: land.scapes[scape_num].raster[int(ind.y), int(ind.x)] for i, ind in self.individs.items()}
        else:
            if scape_num is None:
                hab = {i: [sc.raster[int(ind.y), int(ind.x)] for sc in land.scapes.values()] for i, ind in
                       self.individs.items() if i in individs}
            else:
                hab = {i: land.scapes[scape_num].raster[int(ind.y), int(ind.x)] for i, ind in self.individs.items() if
                       i in individs}

        for i, v in self.individs.items():
            v.habitat = hab[i]

    # method to get all individs' habitat values
    def get_habitat(self):
        return ({i: ind.habitat for i, ind in self.individs.items()})

    # is both args set to None, does same as get_habitat
    def get_habitat_by_land_ind(self, scape_num=None, individs=None):
        if individs is None:
            if scape_num is None:
                return ({i: ind.habitat for i, ind in self.individs.items()})
            else:
                return ({i: ind.habitat[scape_num] for i, ind in self.individs.items()})
        else:
            if scape_num is None:
                return ({i: ind.habitat for i, ind in self.individs.items() if i in individs})
            else:
                return ({i: ind.habitat[scape_num] for i, ind in self.individs.items() if i in individs})

    def get_age(self, individs=None):
        if individs != None:
            return {k: ind.age for k, ind in self.individs.items() if k in individs}
        else:
            return {k: ind.age for k, ind in self.individs.items()}

    def get_genotype(self, locus, return_format='mean', individs=None, by_dominance=False):

        if individs == None:
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

    # convenience method for calling selection.get_phenotype() on this pop
    def get_phenotype(self, trait, individs=None):
        if individs != None:
            return ({i: v for i, v in selection.get_phenotype(self, trait).items() if i in individs})
        else:
            return (selection.get_phenotype(self, trait))

    def get_fitness(self):
        return selection.get_fitness(self)

    def hist_fitness(self):
        plt.hist(selection.get_fitness(self).values())

    def get_single_trait_fitness(self, trait):
        return selection.get_single_trait_fitness(self, trait)

    def get_dom(self, locus):
        return {locus: self.genomic_arch.h[locus]}

    def get_coords(self, individs=None):
        if individs != None:
            return ({k: (float(ind.x), float(ind.y)) for k, ind in self.individs.items() if k in individs})
        else:
            return ({k: (float(ind.x), float(ind.y)) for k, ind in self.individs.items()})

    def get_x_coords(self, individs=None):
        if individs != None:
            return ({k: (float(ind.x)) for k, ind in self.individs.items() if k in individs})
        else:
            return ({k: (float(ind.x)) for k, ind in self.individs.items()})

    def get_y_coords(self, individs=None):
        if individs != None:
            return ({k: (float(ind.y)) for k, ind in self.individs.items() if k in individs})
        else:
            return ({k: (float(ind.y)) for k, ind in self.individs.items()})

    def show(self, land, scape_num=None, color='black', colorbar=True, markersize=25, im_interp_method='nearest',
             alpha=False):
        # if land != None:
        if scape_num != None:
            land.scapes[scape_num].show(colorbar=colorbar, im_interp_method=im_interp_method, pop=True)
        else:
            land.show(colorbar=colorbar, im_interp_method=im_interp_method, pop=True)

        c = np.array(self.get_coords().values())
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
        c = np.array(self.get_coords(individs).values())
        # NOTE: subtract 0.5 to line up points with imshow grid; see note in the pop.show() definition for details
        x = c[:, 0] - 0.5
        y = c[:, 1] - 0.5
        plt.scatter(x, y, s=markersize, c=color, alpha=alpha);
        plt.xlim(-0.6, land.dims[1] - 0.4);
        plt.ylim(-0.6, land.dims[0] - 0.4)

    # NOTE: perhaps worth figuring out how to label with the individual number!!
    # for k, coord_pair in coords.items():
    # ax = mpl.pyplot.plot(coord_pair[0], coord_pair[1], 'ko', scalex = False, scaley = False, color = color, markersize = 8.5)

    def show_density(self, land, window_width=None, normalize_by='census', max_1=False, color='black', markersize=40,
                     alpha=0.5):
        dens = self.calc_density(land, window_width=window_width, normalize_by=normalize_by, max_1=max_1)
        dens.show(im_interp_method='nearest', pop=True)

        c = np.array(self.get_coords().values())
        # NOTE: subtract 0.5 to line up points with imshow grid; see note in the pop.show() definition for details
        x = c[:, 0] - 0.5
        y = c[:, 1] - 0.5
        plt.scatter(x, y, s=markersize, c=color, alpha=alpha);
        plt.xlim(-0.6, land.dims[1] - 0.4);
        plt.ylim(-0.6, land.dims[0] - 0.4)

    # ax = mpl.pyplot.plot([i[0] - 0.5 for i in c.values()], [i[1] - 0.5 for i in c.values()], 'ko', scalex = False, scaley = False, color = color, markersize = 8.5)

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
                c = np.array(self.get_coords(inds).values())
                x = c[:, 0] - 0.5
                y = c[:, 1] - 0.5
                plt.scatter(x, y, s=markersize, c=colors[n], alpha=alpha);
                plt.xlim(-0.6, land.dims[1] - 0.4);
                plt.ylim(-0.6, land.dims[0] - 0.4)

    # method for plotting individuals colored by their phenotypes for a given trait
    def show_phenotype(self, trait, land, scape_num=None, im_interp_method='nearest', markersize=65, alpha=1):

        from collections import OrderedDict as OD

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
        z = self.get_phenotype(trait)

        data = OD({i: (c[i][0] - 0.5, c[i][1] - 0.5, z[i]) for i in self.individs.keys()}).values()

        plt.scatter([i[0] for i in data], [i[1] for i in data], s=markersize, c=[i[2] for i in data], cmap=cmap,
                    alpha=alpha)
        plt.xlim(-0.6, land.dims[1] - 0.4)
        plt.ylim(-0.6, land.dims[0] - 0.4)
        # plt.scatter(x,y, s = markersize, c = [z[i] for i in inds], cmap = cmap, alpha = alpha);plt.xlim(-0.6, land.dims[1]-0.4); plt.ylim(-0.6, land.dims[0]-0.4)

    # method for plotting individuals colored and sized by their overall fitnesses
    def show_fitness(self, land, scape_num=None, im_interp_method='nearest', min_markersize=60, alpha=1):

        from collections import OrderedDict as OD

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
        min_fit = np.product([1 - np.atleast_2d(t.phi).min() for t in self.genomic_arch.traits.values()])

        # generate an evenly spaced range of the possible fitness values
        fit_vals = np.linspace(min_fit, 1, 50)

        # use index of closest possible fitness val to get a markersize differential (to be added to min markersize) for each individual
        markersize_differentials = {i: 3 * np.abs(fit_vals - w[i]).argmin() for i in self.individs.keys()}

        data = OD({i: (c[i][0] - 0.5, c[i][1] - 0.5, w[i], markersize_differentials[i]) for i in
                   self.individs.keys()}).values()

        plt.scatter([i[0] for i in data], [i[1] for i in data], s=[min_markersize + i[3] for i in data],
                    c=[i[2] for i in data], cmap=cmap, alpha=alpha)
        plt.xlim(-0.6, land.dims[1] - 0.4)
        plt.ylim(-0.6, land.dims[0] - 0.4)

        # plt.scatter(x,y, s = [min_markersize + diff for diff in markersize_differentials] , c = [fits[i] for i in inds], cmap = cmap, alpha = alpha);plt.xlim(-0.6, land.dims[1]-0.4); plt.ylim(-0.6, land.dims[0]-0.4)

    # method for plotting individuals colored by their phenotypes for a given trait, sized by their fitness
    def show_single_trait_fitness(self, trait, land, scape_num=None, im_interp_method='nearest', min_markersize=60,
                                  alpha=1):

        from collections import OrderedDict as OD

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
        markersize_differentials = {i: 3 * np.abs(fit_vals - w[i]).argmin() for i in self.individs.keys()}

        z = self.get_phenotype(trait)

        data = OD({i: (c[i][0] - 0.5, c[i][1] - 0.5, w[i], z[i], markersize_differentials[i]) for i in
                   self.individs.keys()}).values()

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
        plt.hist([ind.age for ind in self.individs.values() if ind.sex == 0], orientation='horizontal', color='pink',
                 alpha=0.6)
        plt.hist([ind.age for ind in self.individs.values() if ind.sex == 1], orientation='horizontal', color='skyblue',
                 alpha=0.6)

    def show_pop_growth(self, params):
        T = range(len(self.Nt))
        x0 = params['N'] / self.K.raster.sum()
        plt.plot(T, [demography.logistic_soln(x0, params['R'], t) * self.K.raster.sum() for t in T], color='red')
        plt.plot(T, self.Nt, color='blue')
        plt.xlabel('t')
        plt.ylabel('N(t)')

    # method to plot (or add to an open plot) individuals' IDs
    def plot_ind_ids(self):
        [plt.text(v[0] - 0.5, v[1] - 0.5, i) for i, v in self.get_coords().items()]

    def pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)


# --------------------------------------
# FUNCTIONS ---------------------------
# --------------------------------------


def create_population(genomic_arch, land, params, burn=False):
    # grab necessary params from params dict

    N = params['N']

    dims = params['dims']

    size = params['size']

    T = params['T']

    assert dims.__class__.__name__ in ['tuple', 'list'], "dims should be expressed as a tuple or a list"
    individs = dict()
    for i in range(N):
        # use individual.create_individual to simulate individuals and add them to the population
        ind = individual.create_individual(genomic_arch, dims)
        individs[i] = ind
        land.mating_grid.add(ind)

    pop = Population(N=N, individs=individs, genomic_arch=genomic_arch, size=size, T=T)

    # get initial habitat values
    pop.set_habitat(land)

    return pop


# function for reading in a pickled pop
def load_pickled_pop(filename):
    import cPickle
    with open(filename, 'rb') as f:
        pop = cPickle.load(f)

    return pop
