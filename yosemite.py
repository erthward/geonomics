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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
import os
import time

# set some image params
img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
           'img/final/')
ax_fontdict = {'fontsize': 12,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': 15,
                'name': 'Bitstream Vera Sans'}


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
mid_row_height = 0.1
scnd_col_width = 0.1
fig = plt.figure(figsize=(6.75 + scnd_col_width * 6.75,
                          5.4 + mid_row_height * 5.4))  # constrained_layout=True)
gs = fig.add_gridspec(5, 6,
                      width_ratios=[1, scnd_col_width, 1, 1, 1, 1],
                      height_ratios=[1, 1, mid_row_height, 1, 1])


# define a little function for adding text to each map in the same spot
def add_text_label(text, x=-120.15, y=38.3, size=12, color='black', lw=4,
                   fg='w'):
    txt = plt.text(x, y, text, size=size, color=color)
    txt.set_path_effects([patheffects.withStroke(linewidth=lw, foreground=fg)])


# burn in, then plot starting population, on both rasters
mod.walk(20000, 'burn')

# plot the two rasters before climate change
# axes for tmp before climate change
ax1 = fig.add_subplot(gs[0, 0])
mod.plot(lyr=0, ticks=False, cbar='force')
add_text_label('A')

# cbar = plt.colorbar()
# cbar.set_label('temperature', size=10, rotation=270)

# axes for hab before climate change
ax2 = fig.add_subplot(gs[1, 0])
mod.plot(lyr=1, ticks=False, cbar='force')
add_text_label('B')

# walk for 500 timesteps, then plot population and neighborhood-meaned
# raster before climate change starts
mod.walk(500)

# axes for phenotype plot before climate change
ax3 = fig.add_subplot(gs[:2, 2:4])
mod.plot(lyr=0, cbar=False)
coords = mod.comm[0]._get_plot_coords()
ax3.scatter(x=coords[:, 0], y=coords[:, 1],
            c=[i.z[0] for i in mod.comm[0].values()],
            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
            alpha=1, vmin=0, vmax=1)
add_text_label('C')

# axes for neighborhood mean rast before climate change
ax4 = fig.add_subplot(gs[:2, 4:])
ax4.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
ax4.set_xticks([])
ax4.set_yticks([])
add_text_label('D', 10, 10)

# walk for 100 more timesteps, then plot again,
# at end of climate-change period
mod.walk(100)

# end timer
stop = time.time()
tot_time = stop-start

# axes for tmp after climate change
ax5 = fig.add_subplot(gs[3, 0])
mod.plot(lyr=0, ticks=False, cbar='force')
add_text_label('E')

# axes for hab after climate change
ax6 = fig.add_subplot(gs[4, 0])
mod.plot(lyr=1, ticks=False, cbar='force')
add_text_label('F')

# axes for phenotype plot after climate change
ax7 = fig.add_subplot(gs[3:, 2:4])
mod.plot(lyr=0, cbar=False)
coords = mod.comm[0]._get_plot_coords()
plt.scatter(x=coords[:, 0], y=coords[:, 1],
            c=[i.z[0] for i in mod.comm[0].values()],
            s=ms, cmap=z_cmap, linewidth=0.5, edgecolor='black',
            alpha=1, vmin=0, vmax=1)
add_text_label('G')

# axes for neighborhood mean rast after climate change
ax8 = fig.add_subplot(gs[3:, 4:])
ax8.imshow(calc_neighborhood_mean_phenotype(mod), cmap=z_cmap)
ax8.set_xticks([])
ax8.set_yticks([])
add_text_label('H', 10, 10)

# tweak figure layout
# fig.tight_layout()
# plt.subplots_adjust(left=0,
#                    bottom=0.06,
#                    right=1,
#                    top=0.98,
#                    wspace=0,
#                    hspace=0.27)

# save plot
plt.show()
plt.savefig(os.path.join(img_dir, 'YOSEMITE_time_series.pdf'),
            format='pdf', dpi=1000)

# print out time
print("\n\nModel ran in %0.2f seconds." % tot_time)
