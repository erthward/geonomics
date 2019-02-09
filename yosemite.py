#!/usr/bin/python
# yosemite_example.py

import geonomics as gnx
import utils.spatial as spt

import numpy as np
import utils.io as io
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# a little fn to calculate a future-temperate raster from the 30-year normals
# raster
def calc_fut_tmp(tmp, increment=2):
    fut_tmp = abs(tmp - tmp.max())
    fut_tmp = fut_tmp/fut_tmp.max()
    fut_tmp = tmp + increment + (increment * fut_tmp)
    return fut_tmp


# little fn to calculate a habitat raster that's 1 at the center of env-var's
# range and drops toward 0 at extremes
def calc_hab_rast(tmp, l_lim=7, u_lim=11):
    hab = np.ones(tmp.shape)
    hab[tmp < l_lim] = (1 - (tmp[tmp < l_lim] - l_lim) / (tmp.min() - l_lim))
    hab[tmp > u_lim] = (1 - (u_lim - tmp[tmp > u_lim]) / (u_lim - tmp.max()))
    # hab = 1-abs(env_rast-0.5)
    # hab = (hab - hab.min()) / (hab.max() - hab.min())
    return(hab)


def calc_neighborhood_mean_phenotype(mod, window_width=8):
    # array to store each cell's mean phenotype value
    mean_z = np.ones(mod.land.dim)
    # calc half-window_width
    hww = int(window_width / 2)
    # loop over cells
    for i in range(mod.land.dim[1]):
        for j in range(mod.land.dim[0]):
            # get window around each cell
            i_min = max(i - hww, 0)
            i_max = min(i + hww, mod.land.dim[1])
            j_min = max(j - hww, 0)
            j_max = min(j + hww, mod.land.dim[0])
            # get all phenotypes in the window
            zs_in_window = [i.z for i in mod.comm[0].values(
                ) if ((i_min <= int(i.y) <= i_max) and (
                    j_min <= int(i.x) <= j_max))]
            # average the window's phenotypes and add to mean_z
            # NOTE: if there are no individuals in the window then a NaN is
            # returned
            mean_z[i, j] = np.mean(zs_in_window)

    return mean_z


# create colormap to match the phenotype colors
colors = ['#3C22B4', '#80A6FF', '#FFFFFF']
z_cmap = LinearSegmentedColormap.from_list('my_cmap', colors, N=50)
z_cmap.set_bad(color='black')

print('\nPREPARING MODEL...\n\n')

# read in the params
params = gnx.read_parameters_file('./example/yosemite/yosemite_params.py')

# manually get the temperature raster
tmp = io._read_raster('./example/yosemite/yosemite_30yr_normals_90x90.tif')[0]

# calculate a future-change raster, where change is fastest at highest
# elevations (2 Kelvin raster-wide, plus an additional fraction of 2 that is
# greatest at highest elevations)
fut_tmp = abs(tmp - tmp.max())
fut_tmp = fut_tmp/fut_tmp.max()
fut_tmp = tmp + 2 + (2 * fut_tmp)

# set the scale-min and scale-max values for the tmp layer
scale_min = min(tmp.min(), fut_tmp.min())
scale_max = max(tmp.max(), fut_tmp.max())
params.landscape.layers['tmp'].init.file['scale_min_val'] = scale_min
params.landscape.layers['tmp'].init.file['scale_max_val'] = scale_max

# use tmp rasters to create start and end habitat rasters, and add to params
hab = calc_hab_rast(tmp)
fut_hab = calc_hab_rast(fut_tmp)
# and add them to the params
params.landscape.layers['hab'].init.defined['rast'] = hab
params.landscape.layers['hab'].change[0]['change_rast'] = fut_hab

# scale the fut_tmp raster, then add it to the params
fut_tmp, fut_min, fut_max = spt._scale_raster(fut_tmp, scale_min, scale_max)
params.landscape.layers['tmp'].change[0]['change_rast'] = fut_tmp

# create the model
mod = gnx.make_model(params)

# set plotting params
ms = 6

# burn in, then plot starting population, on both rasters
mod.walk(20000, 'burn')
fig = plt.figure()
ax1 = fig.add_subplot(331)
ax1.set_title('starting population\n(genotypes randomly assigned)')
mod.plot_phenotype(0, 0, 0, size=ms)
ax4 = fig.add_subplot(334)
mod.plot_phenotype(0, 0, 1, size=ms)
ax7 = fig.add_subplot(337)
plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
plt.scatter(x=[i.x for i in mod.comm[0].values()],
            y=[i.y for i in mod.comm[0].values()],
            c=[i.z[0] for i in mod.comm[0].values()],
            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
            alpha=1, vmin=0, vmax=1)


# walk for 500 timesteps, then plot again, before climate change starts
mod.walk(500)
ax2 = fig.add_subplot(332)
ax2.set_title('after 500 timesteps,\n(before climate change begins)')
mod.plot_phenotype(0, 0, 0, size=ms)
ax5 = fig.add_subplot(335)
mod.plot_phenotype(0, 0, 1, size=ms)
ax8 = fig.add_subplot(338)
plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
plt.scatter(x=[i.x for i in mod.comm[0].values()],
            y=[i.y for i in mod.comm[0].values()],
            c=[i.z[0] for i in mod.comm[0].values()],
            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
            alpha=1, vmin=0, vmax=1)


# walk for 1000 more timesteps, then plot again,
# at end of climate-change period
mod.walk(1000)
ax3 = fig.add_subplot(333)
ax3.set_title('after 1500 timesteps\n(at end of period of climate change)')
plt.imshow(mod.land[0].rast, cmap='terrain')
cbar = plt.colorbar()
cbar.set_label('environment', rotation=270)
mod.plot_phenotype(0, 0, 0, size=ms)
ax6 = fig.add_subplot(336)
mod.plot_phenotype(0, 0, 1, size=ms)
ax9 = fig.add_subplot(339)
plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
plt.scatter(x=[i.x for i in mod.comm[0].values()],
            y=[i.y for i in mod.comm[0].values()],
            c=[i.z[0] for i in mod.comm[0].values()],
            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
            alpha=1, vmin=0, vmax=1)
cbar = plt.colorbar()
cbar.set_label('phenotype', rotation=270)
fig.suptitle(('Evolutionary response of 100-locus additive trait to climate '
              'change in Yosemite region, N = ~14300 individuals\n'
              'row 1: temperature rasters (i.e. selective environment); '
              'row 2: habitat quality rasters (i.e. carrying capacity);\n'
              'row 3: raster of neighborhood-meaned phenotypes'))
ax1.set_ylabel('temperature rasters')
ax4.set_ylabel('habitat rasters')
ax7.set_ylabel('neighborhood-meaned phenotype')

# TODO:
# add colorbars for phenotype and for rasters
# create some interpolated/rasterized map of mean phenotye

plt.show()

# TODO run analyses?
