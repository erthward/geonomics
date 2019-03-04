#!/usr/bin/python
# yosemite_example.py

import geonomics as gnx

import numpy as np
import utils.io as io
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from matplotlib.colors import LinearSegmentedColormap
import os

# set some image params
img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
           'img/final/')
ax_fontdict = {'fontsize': 12,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': 15,
                'name': 'Bitstream Vera Sans'}

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
#colors = ['#3C22B4', '#80A6FF', '#FFFFFF']
#z_cmap = LinearSegmentedColormap.from_list('my_cmap', colors, N=50)
#z_cmap.set_bad(color='black')
z_cmap = mpl.cm.coolwarm

print('\nPREPARING MODEL...\n\n')

# read in the params
params = gnx.read_parameters_file(('./geonomics/example/yosemite/'
                                   'yosemite_params.py'))

# manually get the temperature raster
tmp = io._read_raster(('./geonomics/example/yosemite/'
                       'yosemite_30yr_normals_90x90.tif'))[0]

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
fut_tmp, fut_min, fut_max = gnx.utils.spatial._scale_raster(fut_tmp,
                                                            scale_min,
                                                            scale_max)
params.landscape.layers['tmp'].change[0]['change_rast'] = fut_tmp

# create the model
mod = gnx.make_model(params)

# set plotting params
ms = 6

# burn in, then plot starting population, on both rasters
mod.walk(20000, 'burn')
fig = plt.figure()
gs = gridspec.GridSpec(4, 3)
gs.update(wspace=0.1, hspace=0.1)
ax1 = plt.subplot(gs[0, 0])
#mod.plot_phenotype(0, 0, 0, size=ms)
plt.imshow(mod.land[0].rast, cmap='coolwarm')
plt.scatter(x=[i.x for i in mod.comm[0].values()],
            y=[i.y for i in mod.comm[0].values()],
            c=[i.z[0] for i in mod.comm[0].values()],
            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
            alpha=1, vmin=0, vmax=1)
ax4 = plt.subplot(gs[1, 0])
#mod.plot_phenotype(0, 0, 1, size=ms)
plt.imshow(mod.land[0].rast, cmap='coolwarm')
ax7 = plt.subplot(gs[2, 0])
plt.imshow(mod.land[1].rast, cmap='BrBG_r')
ax10 = plt.subplot(gs[3, 0])
plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
#plt.scatter(x=[i.x for i in mod.comm[0].values()],
#            y=[i.y for i in mod.comm[0].values()],
#            c=[i.z[0] for i in mod.comm[0].values()],
#            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
#            alpha=1, vmin=0, vmax=1)


# walk for 500 timesteps, then plot again, before climate change starts
mod.walk(500)
ax2 = plt.subplot(gs[0, 1])
plt.imshow(mod.land[0].rast, cmap='coolwarm')
plt.scatter(x=[i.x for i in mod.comm[0].values()],
            y=[i.y for i in mod.comm[0].values()],
            c=[i.z[0] for i in mod.comm[0].values()],
            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
            alpha=1, vmin=0, vmax=1)
ax5 = plt.subplot(gs[1, 1])
#mod.plot_phenotype(0, 0, 1, size=ms)
plt.imshow(mod.land[0].rast, cmap='coolwarm')
ax8 = plt.subplot(gs[2, 1])
plt.imshow(mod.land[1].rast, cmap='BrBG_r')
ax11 = plt.subplot(gs[3, 1])
plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)

#mod.plot_phenotype(0, 0, 0, size=ms)
#ax5 = plt.subplot(gs[1, 1])
#mod.plot_phenotype(0, 0, 1, size=ms)
#ax8 = plt.subplot(gs[2, 1])
#plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
#plt.scatter(x=[i.x for i in mod.comm[0].values()],
#            y=[i.y for i in mod.comm[0].values()],
#            c=[i.z[0] for i in mod.comm[0].values()],
#            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
#            alpha=1, vmin=0, vmax=1)


# walk for 1000 more timesteps, then plot again,
# at end of climate-change period
mod.walk(1000)
ax3 = plt.subplot(gs[0, 2])
plt.imshow(mod.land[0].rast, cmap='coolwarm')
plt.scatter(x=[i.x for i in mod.comm[0].values()],
            y=[i.y for i in mod.comm[0].values()],
            c=[i.z[0] for i in mod.comm[0].values()],
            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
            alpha=1, vmin=0, vmax=1)
cbar = plt.colorbar()
cbar.set_label('phenotype', rotation=270)
ax6 = plt.subplot(gs[1, 2])
plt.imshow(mod.land[0].rast, cmap='coolwarm')
cbar = plt.colorbar()
cbar.set_label('temperature', rotation=270)
#mod.plot_phenotype(0, 0, 1, size=ms)
ax9 = plt.subplot(gs[2, 2])
plt.imshow(mod.land[1].rast, cmap='BrBG_r')
cbar = plt.colorbar()
cbar.set_label('habitat suitability', rotation=270)
ax12 = plt.subplot(gs[3, 2])
plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
cbar = plt.colorbar()
cbar.set_label('phenotype', rotation=270)


#ax3 = plt.subplot(gs[0, 2])
#plt.imshow(mod.land[0].rast, cmap='terrain')
#cbar = plt.colorbar()
#cbar.set_label('environment', rotation=270)
#mod.plot_phenotype(0, 0, 0, size=ms)
#ax6 = plt.subplot(gs[1, 2])
#mod.plot_phenotype(0, 0, 1, size=ms)
#ax9 = plt.subplot(gs[2, 2])
#plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
#plt.scatter(x=[i.x for i in mod.comm[0].values()],
#            y=[i.y for i in mod.comm[0].values()],
#            c=[i.z[0] for i in mod.comm[0].values()],
#            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
#            alpha=1, vmin=0, vmax=1)
# fig.suptitle(('Evolutionary response of 100-locus additive trait to climate '
#              'change in Yosemite region, N = ~14300 individuals\n'
#              'row 1: temperature rasters (i.e. selective environment); '
#              'row 2: habitat quality rasters (i.e. carrying capacity);\n'
#              'row 3: raster of neighborhood-meaned phenotypes'))

ax1.set_title('starting population\n(genotypes randomly assigned)')
ax2.set_title('after 500 timesteps,\n(before climate change begins)')
ax3.set_title('after 1500 timesteps\n(after climate change)')
ax1.set_ylabel(('population (plotted on temperature rasters,\ncolored by '
                'phenotype)'), fontdict=ax_fontdict)
ax4.set_ylabel('temperature rasters', fontdict=ax_fontdict)
ax7.set_ylabel('habitat rasters', fontdict=ax_fontdict)
ax10.set_ylabel('neighborhood-meaned phenotype', fontdict=ax_fontdict)

plt.show()
plt.savefig(os.path.join(img_dir, 'YOSEMITE_time_series.png'))
# TODO run analyses?
