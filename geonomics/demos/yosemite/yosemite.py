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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import time

# set the directory where the gif images should be stored
gif_dir = '/home/drew/Desktop/stuff/berk/research/projects/sim/gif_dir/'

# set some image params
img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
           'img/final/')
ax_fontdict = {'fontsize': 12,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': 15,
                'name': 'Bitstream Vera Sans'}

# set the amount of time before and after climate change
t_before_cc = 500
t_after_cc = 100


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
params = gnx.read_parameters_file(('./geonomics/demos/yosemite/'
                                   'yosemite_params.py'))

# create the model
mod = gnx.make_model(params, verbose=True)

# set plotting params
ms = 6

# set up the multipanel plot with gridspec, for the methods-paper figure
mid_row_height = 0.1
scnd_col_width = 0.1
fig = plt.figure(figsize=(6.75 + scnd_col_width * 6.75,
                          5.4 + mid_row_height * 5.4))
#                         constrained_layout=True)
gs = fig.add_gridspec(5, 6,
                      width_ratios=[1, scnd_col_width, 1, 1, 1, 1],
                      height_ratios=[1, 1, mid_row_height, 1, 1])


# define a little function for adding text to each map in the same spot
def add_text_label(text, x=-120.15, y=38.3, size=12, color='black', lw=4,
                   fg='w'):
    txt = plt.text(x, y, text, size=size, color=color)
    txt.set_path_effects([patheffects.withStroke(linewidth=lw, foreground=fg)])


# function for creating a saving images that imagemagic will stitch into a GIF
def save_gif_img(t):
    # set up second figure, for a GIF animation
    mid_col_width = 0.8
    fig2 = plt.figure(figsize=(6.75 + mid_col_width * 6.75, 5.4))
    ax2_gs = fig.add_gridspec(1, 3, width_ratios=[1, scnd_col_width, 1])
    # fig2, axs2 = plt.subplots(1, 2)
    # fig2.set_tight_layout(True)
    ax2_1 = fig2.add_subplot(ax2_gs[0, 0])
    ax2_2 = fig2.add_subplot(ax2_gs[0, 2])
    fig2.suptitle('Time step: %i' % (t+1), size=20)
    ax2_1.set_title('Temperature', size=12)
    ax2_2.set_title('Habitat suitability', size=12)
    ax2_1.set_xlabel('lon', size=8)
    ax2_1.set_ylabel('lat', size=8)
    ax2_2.set_xlabel('lon', size=8)
    ax2_2.set_ylabel('lat', size=8)
    ax2_1.set_aspect(1)
    ax2_2.set_aspect(1)
    # plot the population on both the temp and hab rasters
    land1 = ax2_1.pcolormesh(mod.land._x_cell_bds,
                             mod.land._y_cell_bds,
                             mod.land[0]._get_rast_in_native_units(),
                             cmap='RdBu_r',
                             vmin=mod.land[0]._scale_min,
                             vmax=mod.land[0]._scale_max)
    xticks, xticklabs, yticks, yticklabs = mod.land[0]._get_coord_ticks()
    ax2_1.set_xticklabels(labels=xticklabs, size=6)
    ax2_1.set_yticklabels(labels=yticklabs, size=6)
    divider = make_axes_locatable(ax2_1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar_ticks = mod.land[0]._get_cbar_ticks_and_minmax_scaled_vals()
    cbar = plt.colorbar(land1, cax=cax, ticks=np.linspace(cbar_ticks[1],
                                                          cbar_ticks[2], 5))
    cbar.ax.set_yticklabels(cbar_ticks[0], size=6)
    cbar.set_label(mod.land[0].units, rotation=270, labelpad=5, y=0.5,
                   size=7)
    coords = mod.comm[0]._get_plot_coords()
    xs = coords[:, 0]
    ys = coords[:, 1]
    zs = mod.comm[0]._get_z()[:, 0]
    ax2_1.scatter(xs, ys, c=zs, s=1, cmap='RdBu_r')
    land2 = ax2_2.pcolormesh(mod.land._x_cell_bds,
                             mod.land._y_cell_bds,
                             mod.land[1]._get_rast_in_native_units(),
                             cmap='BrBG_r',
                             vmin=mod.land[1]._scale_min,
                             vmax=mod.land[1]._scale_max)
    xticks, xticklabs, yticks, yticklabs = mod.land[1]._get_coord_ticks()
    ax2_2.set_xticklabels(labels=xticklabs, size=6)
    ax2_2.set_yticklabels(labels=yticklabs, size=6)
    divider = make_axes_locatable(ax2_2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar_ticks = mod.land[1]._get_cbar_ticks_and_minmax_scaled_vals()
    cbar = plt.colorbar(land2, cax=cax, ticks=np.linspace(cbar_ticks[1],
                                                          cbar_ticks[2], 5))
    cbar.ax.set_yticklabels(cbar_ticks[0], size=6)
    cbar.set_label(mod.land[1].units, rotation=270, labelpad=15, y=0.5,
                   size=10)
    ax2_2.scatter(xs, ys, c='black', s=0.5)

    # save the figure
    fig2.savefig(os.path.join(gif_dir, 'gif_img_%i.jpg' % (t + 1)))

    # manually close the figure
    plt.close(fig2)


# burn in, then plot starting population, on both rasters
mod.walk(20000, 'burn')

# run for the first 500 timsteps, before climate change
for _ in range(t_before_cc):
    mod.walk(1)
    save_gif_img(mod.t)

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
neigh_mean = calc_neighborhood_mean_phenotype(mod)
ax4.imshow(neigh_mean, cmap=z_cmap)
# add a contour at phenotype = 0.5
ax4.contour(neigh_mean, colors=['#17161a'], alpha=0.5,
            levels=np.array([0.5]))
ax4.set_xticks([])
ax4.set_yticks([])
add_text_label('D', 10, 10)

# walk for 100 more timesteps, then plot again,
# at end of climate-change period
for _ in range(t_after_cc):
    mod.walk(1)
    save_gif_img(mod.t)

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
neigh_mean = calc_neighborhood_mean_phenotype(mod)
ax8.imshow(neigh_mean, cmap=z_cmap)
ax4.contour(neigh_mean, colors=['#17161a'], alpha=0.5,
            levels=np.array([0.5]))
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

# create and save a population-size plot
mod.walk(50)
fig3 = plt.figure()
ax3 = fig.add_subplot(111)
burn_len = mod.burn_t
line_height = int(10000*np.ceil(max(mod.comm[0].Nt)/10000))
plt.plot(range(len(mod.comm[0].Nt)), mod.comm[0].Nt)
plt.plot([burn_len, burn_len], [0, line_height], c='red')
plt.plot([burn_len+500, burn_len+500], [0, line_height], c='red')
plt.plot([burn_len+600, burn_len+600], [0, line_height], c='red')
chng_yrs = [int(x[:3]) for x in os.listdir(('./geonomics/examples/yosemite/'
                                            'yosemite_lyrs/ppt'))]
for yr in chng_yrs:
    plt.plot([burn_len+yr, burn_len+yr], [0, line_height], ':r', linewidth=0.5)
ax3.set_label('time (time steps/years)')
ax3.set_ylabel('total population size (individuals)')
plt.savefig(os.path.join(img_dir, 'YOSEMITE_time_series_POP_SIZE.png'),
            format='png')

# create the GIF using imagemagick
os.system('cd %s' % gif_dir)
os.system(('cd %s; convert -delay 5 -loop 0 `ls -v` '
           '../yosemite.gif; cd ..') % gif_dir)
os.system('..')

# print out time
print("\n\nModel ran in %0.2f seconds." % tot_time)
