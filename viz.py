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


def show_rasters(land, scape_num = None, colorbar = True, im_interp_method = 'nearest', cmap = 'terrain', plt_lims = None, mask_val = None):
    #if a figure is already open, force colorbar to False
    if plt.get_fignums():
        colorbar = False

    #if just a numpy.ndarray or a Landscape (not a Landscape_Stack) is provided, or if just a single raster is desired, grab the raster into a list
    if str(type(land)) == "<class 'numpy.ndarray'>":
        rasters = [land]
    elif str(type(land)) == "<class 'landscape.Landscape'>":
        rasters = [land.raster]
    elif str(type(land)) == "<class 'landscape.Landscape_Stack'>":
        if scape_num is not None:
            rasters = [land.scapes[scape_num].raster]
        #else just create a list of all rasters
        else:
            rasters = [land.scapes[i].raster for i in range(len(land.scapes))]

    if type(cmap) == str:
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
        #if zoom is not None:
        #    min_i, max_i = zoom[0]
        #    min_j, max_j = zoom[1]
        #    rasters[n] = np.array([row[min_j:max_j] for row in rasters[n][min_i:max_i]])
        plt.imshow(rasters[n], interpolation=im_interp_method, cmap=cmaps[n], vmin=vmin, alpha = alphas[n])
        if plt_lims is not None:
            plt.xlim(plt_lims[0])
            plt.ylim(plt_lims[1])
        #and their colorbars, if requested (but for only the first two rasters maximum, 
        #since the second and onward share the same palette)
        if colorbar and n < 2:
            clb = plt.colorbar(boundaries=np.linspace(0, max(rasters[n].max(), 1), 51))
            clb.ax.set_title('Scape: %i' % n)
            ax = clb.ax
            text = ax.title
            font = mpl.font_manager.FontProperties(family='arial', style='normal', size=10)
            text.set_font_properties(font)


def show_points(points, scape_num=None, color='black', edge_color='face', text_color='purple', linewidth=0.5,
        pt_cmap=None, size=25, text_size=12, alpha=False, text=None, plt_lims=None):
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
    if pt_cmap == 'terrain': 
        from matplotlib.colors import LinearSegmentedColormap
        colors = ['#3C22B4', '#80A6FF', '#FFFFFF']
        # colors to match matplotlib colormap 'terrain' palette extremes, but with hybrid a mix of the extremes
        # rather than the yellow at the middle of the palette, for nicer viewing
        cmap = LinearSegmentedColormap.from_list('my_cmap', colors, N=50)
    
        plt.scatter(x, y, s=size, c=color, cmap=cmap, linewidth=linewidth, edgecolor=edge_color, alpha=alpha)
    else:
        plt.scatter(x, y, s=size, c=color, linewidth=linewidth, edgecolor=edge_color, alpha=alpha);

    #add text, if requested
    if text is not None:
        show_text = []
        for n,t in enumerate(text):
            if plt_lims is not None:
                if plt_lims[0][0] <= x[n] <= plt_lims[0][1] and plt_lims[1][0] <= y[n] <= plt_lims[1][1]:
                    show_text.append((x[n], y[n], t))
            else:
                show_text.append((x[n], y[n], t))
        [plt.text(*item, color=text_color, size=text_size, alpha=alpha) for item in show_text];

    if plt_lims is not None and len(plt_lims) == 2 and [len(item) for item in plt_lims] == [2,2]:
        plt.xlim(plt_lims[0])
        plt.ylim(plt_lims[1])
    else:
        print("plt_lims appears not to be a valid argument (i.e. a 2-tuple of 2-tuples)")


def get_scape_plt_lims(land):
    xlim = (-1, land.dims[1])
    ylim = (-1, land.dims[0])
    lims = (xlim, ylim)
    return(lims)

def get_zoom_plt_lims(x, y, zoom_width):
    #get zoom-half-width
    zhw = zoom_width/2
    xlim = (x- zhw, x+zhw)
    ylim = (y- zhw, y+zhw)
    lims = (xlim, ylim)
    return(lims)
   
def get_plt_lims(land=None, x=None, y=None, zoom_width=None):
    if zoom_width is not None and x is not None and y is not None:
        plt_lims = get_zoom_plt_lims(x, y, zoom_width) 
    else: 
        plt_lims = get_scape_plt_lims(land)
    return(plt_lims)







#### FROM POPULATION.PY

# method for plotting individuals colored by their phenotypes for a given trait
def show_phenotype(self, trait, land, scape_num=None, im_interp_method='nearest', marker_size=65, alpha=1):

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

    plt.scatter([i[0] for i in data], [i[1] for i in data], s=marker_size, c=[i[2] for i in data], cmap=cmap,
                linewidth=1, edgecolor='black',
                alpha=alpha)
    plt.xlim(-0.6, land.dims[1] - 0.4)
    plt.ylim(-0.6, land.dims[0] - 0.4)
    # plt.scatter(x,y, s = marker_size, c = [z[i] for i in i.items()}nds], cmap = cmap, alpha = alpha);plt.xlim(-0.6, land.dims[1]-0.4); plt.ylim(-0.6, land.dims[0]-0.4)

# method for plotting individuals colored and sized by their overall fitnesses
def show_fitness(self, land, scape_num=None, im_interp_method='nearest', min_marker_size=60, alpha=1):

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

    # use index of closest possible fitness val to get a marker_size differential (to be added to min marker_size) for each individual
    marker_size_differentials = {i: 3 * np.abs(fit_vals - w[n]).argmin() for n,i in enumerate(self.individs.keys())}

    data = list(OD({i: (c[i][0] - 0.5, c[i][1] - 0.5, w[i], marker_size_differentials[i]) for i in
                    self.individs.keys()}).values())

    plt.scatter([i[0] for i in data], [i[1] for i in data], s=[min_marker_size + i[3] for i in data],
                c=[i[2] for i in data], cmap=cmap, alpha=alpha)
    plt.xlim(-0.6, land.dims[1] - 0.4)
    plt.ylim(-0.6, land.dims[0] - 0.4)


# method for plotting individuals colored by their phenotypes for a given trait, sized by their fitness
def show_single_trait_fitness(self, trait, land, scape_num=None, im_interp_method='nearest', min_marker_size=60,
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

    # use index of closest possible fitness val to get a marker_size differential (to be added to min marker_size) for each individual
    marker_size_differentials = [3 * np.abs(fit_vals - w[i]).argmin() for i in range(self.census())]

    z = [v[trait] for v in self.get_phenotype()]

    data = list(OD({i: (c[n][0] - 0.5, c[n][1] - 0.5, w[n], z[n], marker_size_differentials[n]) for n,i in
                    enumerate(self.individs.keys())}).values())

    # separate colormap to color marker edges from black (fit = 1) to white (fit = 0) through red
    inside_marker_cmap = mpl.cm.get_cmap('RdYlGn')

    plt.scatter([i[0] for i in data], [i[1] for i in data], s=[min_marker_size + i[4] for i in data],
                c=[i[3] for i in data], cmap=cmap, alpha=alpha)
    plt.scatter([i[0] for i in data], [i[1] for i in data], s=min_marker_size / 3, c=[i[2] for i in data],
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

