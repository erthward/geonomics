#!/usr/bin/python
# sim_sel.py

# flake8: noqa

import geonomics as gnx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import os

# set the directory for saving methods-paper images
img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/'
           'sim/methods_paper/img/final')


# define a function to calculate the mean difference between phenotype
# and environment for a species
def calc_mean_z_e_diff(spp, trait_num=0):
    zs = spp._get_z()[:, trait_num].ravel()
    es = spp._get_e(lyr_num=spp.gen_arch.traits[trait_num].lyr_num)
    mean_diff = np.mean(np.abs(zs - es))
    return mean_diff


# set flag to indicate whether or not I'm timing this demo
# (because if I am then I don't want to be calculating the mean z-e diff
# at each timestep)
not_timing = True


# create data structure to store z-e diffs for both traits
z_e_diffs = {0: [], 1: []}

# start timer
start = time.time()

# RUN THE MODEL WITHOUT LINKAGE
# make the model
mod = gnx.make_model(('./geonomics/demos/simult_select/'
                      'simult_select_params.py'))

# get total runtime from params
T = mod.params.model.T

# run model
mod.walk(T=10000, mode='burn', verbose=True)
for t in range(T):
    if not_timing:
        for trt_num in range(len(mod.comm[0].gen_arch.traits)):
            z_e_diffs[trt_num].append(calc_mean_z_e_diff(mod.comm[0],
                                                         trait_num=trt_num))
    mod.walk(T=1, mode='main', verbose=True)
if not_timing:
    for trt_num in range(len(mod.comm[0].gen_arch.traits)):
        z_e_diffs[trt_num].append(calc_mean_z_e_diff(mod.comm[0],
                                                     trait_num=trt_num))

# stop timer
stop = time.time()
tot_time = stop - start

# plot the resulting species on top of each layer
fig = plt.figure(figsize=(9.25, 4.5))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
im = plt.pcolormesh(np.linspace(0, mod.land.dim[0], mod.land.dim[0]+1),
                    np.linspace(0, mod.land.dim[1], mod.land.dim[1]+1),
                    mod.land[0].rast, cmap='coolwarm')
ax1 = plt.subplot(gs[0])
mod.plot_phenotype(0, 0, 0, size=85)
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('environmental and phenotypic value', rotation=270,
               labelpad=25, y=0.5, size=24)
ax2 = plt.subplot(gs[1])
im = plt.pcolormesh(np.linspace(0, mod.land.dim[0], mod.land.dim[0]+1),
                    np.linspace(0, mod.land.dim[1], mod.land.dim[1]+1),
                    mod.land[1].rast, cmap='BrBG_r')
mod.plot_phenotype(0, 1, 1, size=85)
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('environmental and phenotypic value', rotation=270,
               labelpad=25, y=0.5, size=24)
ax1.set_title('trait 0', size=30)
ax2.set_title('trait 1', size=30)
plt.show()
fig.tight_layout()
# print out time
print("Model ran in %0.2f seconds." % tot_time)

# save the figure
plt.savefig(os.path.join(img_dir, 'sim_sel.pdf'), format='pdf', dpi=1000)


# plot z-e diffs
z_e_fig = plt.figure()
ax = z_e_fig.add_subplot(111)
trt_colors = ['#096075', '#e87d4f']
plt.plot(range(len(z_e_diffs[0])), z_e_diffs[0], trt_colors[0])
plt.plot(range(len(z_e_diffs[1])), z_e_diffs[1], trt_colors[1])
ax.legend(labels=['trait = %i' % trt for trt in range(2)],
          loc='best', fontsize='medium')
ax.set_xlabel('time')
ax.set_ylabel(('mean difference between individuals\' phenotypes and '
               'environmental values'))
plt.show()
