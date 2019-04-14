#!/usr/bin/python
# yosemite_example.py

'''

TODO
1.  Add a "timing" argument, that I can set to false to run the model without
plotting, to get an accurate assessment of run time

2. Debug all the plotting stuff

'''

import geonomics as gnx

import numpy as np
# import utils.io as io
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# from matplotlib.colors import LinearSegmentedColormap
# import os
import time

# set some image params
img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
           'img/final/')
ax_fontdict = {'fontsize': 12,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': 15,
                'name': 'Bitstream Vera Sans'}

# a little fn to calculate a future-temperate raster from the 30-year normals
# raster
# def calc_fut_tmp(tmp, increment=2):
#    fut_tmp = abs(tmp - tmp.max())
#    fut_tmp = fut_tmp/fut_tmp.max()
#    fut_tmp = tmp + increment + (increment * fut_tmp)
#    return fut_tmp
#
#
# little fn to calculate a habitat raster that's 1 at the center of env-var's
# range and drops toward 0 at extremes
# def calc_hab_rast(tmp, l_lim=7, u_lim=11):
#    hab = np.ones(tmp.shape)
#    hab[tmp < l_lim] = (1 - (tmp[tmp < l_lim] - l_lim) / (tmp.min() - l_lim))
#    hab[tmp > u_lim] = (1 - (u_lim - tmp[tmp > u_lim]) / (u_lim - tmp.max()))
#    # hab = 1-abs(env_rast-0.5)
#    # hab = (hab - hab.min()) / (hab.max() - hab.min())
#    return(hab)


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
z_cmap = mpl.cm.RdBu_r

# start timer
start = time.time()

# read in the params
params = gnx.read_parameters_file(('./geonomics/examples/yosemite/'
                                   'yosemite_params.py'))

# create the model
mod = gnx.make_model(params, verbose=True)

# set plotting params
ms = 6

# set up the multipanel plot with gridspec
fig = plt.figure(figsize=(6.75, 7.0875))  # constrained_layout=True)
gs = fig.add_gridspec(4, 5,
                      width_ratios=[1, 1, 1, 1, 1],
                      height_ratios=[1, 1, 1, 1])

# burn in, then plot starting population, on both rasters
mod.walk(20000, 'burn')
# fig = plt.figure()
# gs = gridspec.GridSpec(4, 3)
# gs.update(wspace=0.1, hspace=0.1)
# ax1 = plt.subplot(gs[0, 0])
# mod.plot_phenotype(0, 0, 0, size=ms)
# plt.imshow(mod.land[0].rast, cmap='coolwarm')
# plt.scatter(x=[i.x for i in mod.comm[0].values()],
#            y=[i.y for i in mod.comm[0].values()],
#            c=[i.z[0] for i in mod.comm[0].values()],
#            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
#            alpha=1, vmin=0, vmax=1)
# ax4 = plt.subplot(gs[1, 0])
# mod.plot_phenotype(0, 0, 1, size=ms)
# plt.imshow(mod.land[0].rast, cmap='coolwarm')
# ax7 = plt.subplot(gs[2, 0])
# plt.imshow(mod.land[1].rast, cmap='BrBG_r')
# ax10 = plt.subplot(gs[3, 0])
# plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
# plt.scatter(x=[i.x for i in mod.comm[0].values()],
#            y=[i.y for i in mod.comm[0].values()],
#            c=[i.z[0] for i in mod.comm[0].values()],
#            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
#            alpha=1, vmin=0, vmax=1)

# plot the two rasters before climate change
# axes for tmp before climate change
ax1 = fig.add_subplot(gs[0, 0])
mod.plot(lyr=0)
cbar = plt.colorbar()
cbar.set_label('temperature', size=9, rotation=270)

# axes for hab before climate change
ax2 = fig.add_subplot(gs[1, 0])
mod.plot(lyr=1)
cbar = plt.colorbar()
cbar.set_label('habitat\nsuitability', size=9, rotation=270)

# walk for 500 timesteps, then plot population and neighborhood-meaned
# raster before climate change starts
mod.walk(500)

# axes for phenotype plot before climate change
ax3 = fig.add_subplot(gs[:2, 1:3])
mod.plot(lyr=0)
coords = mod.comm[0]._get_plot_coords()
ax3.scatter(x=coords[:, 0], y=coords[:, 1],
            c=[i.z[0] for i in mod.comm[0].values()],
            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
            alpha=1, vmin=0, vmax=1)

# axes for neighborhood mean rast before climate change
ax4 = fig.add_subplot(gs[:2, 3:])
ax4.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
ax4.set_xticks([])
ax4.set_yticks([])

# mod.plot_phenotype(0, 0, 0, size=ms)
# ax5 = plt.subplot(gs[1, 1])
# mod.plot_phenotype(0, 0, 1, size=ms)
# ax8 = plt.subplot(gs[2, 1])
# plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
# plt.scatter(x=[i.x for i in mod.comm[0].values()],
#            y=[i.y for i in mod.comm[0].values()],
#            c=[i.z[0] for i in mod.comm[0].values()],
#            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
#            alpha=1, vmin=0, vmax=1)


# walk for 100 more timesteps, then plot again,
# at end of climate-change period
mod.walk(100)

# end timer
stop = time.time()
tot_time = stop-start

# axes for tmp after climate change
ax5 = fig.add_subplot(gs[2, 0])
mod.plot(lyr=0)
cbar = plt.colorbar()
cbar.set_label('temperature', size=9, rotation=270)

# axes for hab after climate change
ax6 = fig.add_subplot(gs[3, 0])
mod.plot(lyr=1)
cbar = plt.colorbar()
cbar.set_label('habitat\nsuitability', size=9, rotation=270)

# axes for phenotype plot after climate change
ax7 = fig.add_subplot(gs[2:, 1:3])
mod.plot(lyr=0)
coords = mod.comm[0]._get_plot_coords()
plt.scatter(x=coords[:, 0], y=coords[:, 1],
            c=[i.z[0] for i in mod.comm[0].values()],
            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
            alpha=1, vmin=0, vmax=1)
# cbar = plt.colorbar()
# cbar.set_label('phenotype', rotation=270)
# ax6 = plt.subplot(gs[1, 2])
# plt.imshow(mod.land[0].rast, cmap='coolwarm')
# cbar = plt.colorbar()
# cbar.set_label('temperature', rotation=270)
# mod.plot_phenotype(0, 0, 1, size=ms)
# ax9 = plt.subplot(gs[2, 2])
# plt.imshow(mod.land[1].rast, cmap='BrBG_r')
# cbar = plt.colorbar()
# cbar.set_label('habitat suitability', rotation=270)
# ax12 = plt.subplot(gs[3, 2])

# axes for neighborhood mean rast after climate change
ax8 = fig.add_subplot(gs[2:, 3:])
ax8.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
ax8.set_xticks([])
ax8.set_yticks([])
# cbar = plt.colorbar()
# cbar.set_label('phenotype', rotation=270)


# ax3 = plt.subplot(gs[0, 2])
# plt.imshow(mod.land[0].rast, cmap='terrain')
# cbar = plt.colorbar()
# cbar.set_label('environment', rotation=270)
# mod.plot_phenotype(0, 0, 0, size=ms)
# ax6 = plt.subplot(gs[1, 2])
# mod.plot_phenotype(0, 0, 1, size=ms)
# ax9 = plt.subplot(gs[2, 2])
# plt.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
# plt.scatter(x=[i.x for i in mod.comm[0].values()],
#            y=[i.y for i in mod.comm[0].values()],
#            c=[i.z[0] for i in mod.comm[0].values()],
#            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
#            alpha=1, vmin=0, vmax=1)
# fig.suptitle(('Evolutionary response of 100-locus additive trait to climate '
#              'change in Yosemite region, N = ~14300 individuals\n'
#              'row 1: temperature rasters (i.e. selective environment); '
#              'row 2: habitat quality rasters (i.e. carrying capacity);\n'
#              'row 3: raster of neighborhood-meaned phenotypes'))

# ax1.set_title('starting population\n(genotypes randomly assigned)')
# ax2.set_title('after 500 timesteps,\n(before climate change begins)')
# ax3.set_title('after 1500 timesteps\n(after climate change)')
# ax1.set_ylabel(('population (plotted on temperature rasters,\ncolored by '
#                 'phenotype)'), fontdict=ax_fontdict)
# ax4.set_ylabel('temperature rasters', fontdict=ax_fontdict)
# ax7.set_ylabel('habitat rasters', fontdict=ax_fontdict)
# ax10.set_ylabel('neighborhood-meaned phenotype', fontdict=ax_fontdict)

fig.tight_layout()
plt.show()
# plt.savefig(os.path.join(img_dir, 'YOSEMITE_time_series.pdf'),
#            format='pdf',
#            dpi=1000)

# print out time
print("\n\nModel ran in %0.2f seconds." % tot_time)
