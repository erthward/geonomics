#!/usr/bin/python
# viz.py

'''
##########################################

Module name:          viz


Module contains:
                      - basic visualization functions, to be inherited by classes of model components


Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           06-18-18
Documentation:        URL


##########################################
'''

import numpy as np
import numpy.random as r
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter as C
from operator import itemgetter as ig
from shapely import geometry as g
from operator import itemgetter
from operator import attrgetter
import copy
import sys



def show_rasters(land, scape_num = None, colorbar = True, im_interp_method = 'nearest', cmap = 'terrain', zoom = None, mask_val = None):
    #if a figure is already open, force colorbar to False
    if plt.get_fignums():
        colorbar = False

    #if just a Landscape (not a Landscape_Stack) is provided, or if just a single raster is desired, grab the raster into a list
    if str(type(land)) == "<class 'landscape.Landscape'>": 
        rasters = [land.raster]
    elif str(type(land)) == "<class 'landscape.Landscape_Stack'>":
        if scape_num is not None:
            rasters = [land.scapes[scape_num].raster]
        #else just create a list of all rasters
        else:
            rasters = [land.scapes[i].raster for i in range(len(land.scapes))]

    #get the requested cmap 
    cmap = getattr(plt.cm, cmap)
    #set the minimum plotting value
    vmin = 0
    #mask values below mask_val, if not None
    if mask_val is not None:
        cmap.set_under(color = 'black')
        vmin = mask_val

    #create cmaps and alphas lists, in case multiple rasters are to be plotted
    cmaps = [cmap] + ['bone'] * (len(rasters)-1)
    alphas = [1] + [0.5] * (len(rasters)-1)
    #plot all the rasters...
    for n in range(len(rasters)):
        #pull out the zoomed raster, if requested
        if zoom is not None:
            min_i, max_i = zoom[0]
            min_j, max_j = zoom[1]
            rasters[n] = np.array([row[min_j:max_j] for row in rasters[n][min_i:max_i]])
        plt.imshow(rasters[n], interpolation=im_interp_method, cmap=cmaps[n], vmin=vmin, alpha = alphas[n])
        #and their colorbars, if requested (but for only the first two rasters maximum, 
        #since the second and onward share the same palette)
        if colorbar and n < 2:
            clb = plt.colorbar(boundaries=np.linspace(0, max(rasters[n].max(), 1), 51))
            clb.ax.set_title('Scape: %i' % n)
            ax = clb.ax
            text = ax.title
            font = mpl.font_manager.FontProperties(family='arial', style='normal', size=10)
            text.set_font_properties(font)


def show_points(points, scape_num=None, color='black', terrain_colormap = False, markersize=25, alpha = False, text=None, xlim = None, ylim = None):
    #get the x and y coordinates from the points (and subtract 0.5 to line the points up with the plt.imshow()
    #grid of a landscape raster; imshow plots each pixel centered on its index, but the points then plot on 
    #those indices, so wind up shifted +0.5 on each axis
    x = points[:, 0] - 0.5
    y = points[:, 1] - 0.5
    #handle the alpha value as necessary
    if alpha == True and type(alpha) == bool:
        alpha = 0.6
    elif alpha == True and type(alpha) == float:
        alpha = alpha
    else:
        alpha = 1.0

    #plot the points, as stipulated by arguments
    if text is None:

        if terrain_colormap: 
            from matplotlib.colors import LinearSegmentedColormap
            colors = ['#3C22B4', '#80A6FF', '#FFFFFF']
            # colors to match matplotlib colormap 'terrain' palette extremes, but with hybrid a mix of the extremes
            # rather than the yellow at the middle of the palette, for nicer viewing
            cmap = LinearSegmentedColormap.from_list('my_cmap', colors, N=50)
    
            plt.scatter(x, y, s=markersize, c=color, cmap=cmap, linewidth=1, edgecolor='black', alpha=alpha)
        else:
            plt.scatter(x, y, s=markersize, c=color, alpha=alpha);
    else:
        [plt.text(x[i], y[i], text[i]) for i in text]

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    #plt.ylim(-0.6, land.dims[0] - 0.4)




#### FROM LANDSCAPE.PY
def plot_movement_surf_vectors(self, params, circle=False):
    if params['land']['movement_surf'] == True and self.movement_surf is not None:
        import movement
        movement.plot_movement_surf_vectors(self, self.movement_surf_scape_num, circle=circle)
    else:
        print('\n\nThis Landscape_Stack appears to have no movement surface.\n\n')
        pass

def plot_vonmises_mix_circ_draws(self, i, j, params):
    if self.movement_surf is None:
        print('This landscape stack appears to have no movement surface layer. Function not valid.')
        return
    else:
        scape_num = self.movement_surf_scape_num
        pts = [(np.cos(a), np.sin(a)) for a in [self.movement_surf[i][j]()[0] for n in range(1000)]]
        plt.scatter([pt[0] * 0.5 + i for pt in pts], [pt[1] * 0.5 + j for pt in pts], color='red', alpha=0.2)
        self.scapes[scape_num].zoom(max(i - 10, 0), min(i + 10, self.dims[0]), max(j - 10, 0),
                                    min(j + 10, self.dims[1]))

def plot_vonmises_mix_hist(self, i, j):
    if self.movement_surf is None:
        print('This landscape stack appears to have no movement surface layer. Function not valid.')
        return
    else:
        plt.hist([self.movement_surf[i][j]()[0] for n in range(10000)], bins=100, normed=True, alpha=0.5)

def plot_vonmises_mix_circ_hist(self, i, j, zoom, params, scale_fac=4.5, color='black'):
    if self.movement_surf is None:
        print('This landscape stack appears to have no movement surface layer. Function not valid.')
        return
    else:
        scape_num = land.movement_surf_scape_num
        v, a = np.histogram(r.choice(list(self.movement_surf[i][j]), replace = True, size = 7500), bins=15)
        v = v / float(v.sum())
        a = [(a[n] + a[n + 1]) / 2 for n in range(len(a) - 1)]
        x = [np.cos(a[n]) * 0.5 for n in range(len(a))]
        y = [np.sin(a[n]) * 0.5 for n in range(len(a))]
        x = np.array(x) * v * scale_fac
        y = np.array(y) * v * scale_fac
        [plt.plot((j, (j + x[n])), (i, (i + y[n])), linewidth=2, color=color) for n in range(len(x))]
        self.scapes[scape_num].zoom(max(i - zoom, 0), min(i + zoom, self.dims[0]), max(j - zoom, 0),
                                        min(j + zoom, self.dims[1]), pop=True)







#### FROM POPULATION.PY

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

