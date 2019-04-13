#!/usr/bin/python
# sim_sel.py

import geonomics as gnx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

# start timer
start = time.time()

# RUN THE MODEL WITHOUT LINKAGE
# make the model
mod = gnx.make_model('./geonomics/examples/sim_sel/sim_sel_params.py')
# run it
mod.run(verbose=True)

# stop timer
stop = time.time()
tot_time = stop - start

# plot the resulting species on top of each layer
fig = plt.figure()
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
im = plt.pcolormesh(np.linspace(0, mod.land.dim[0], mod.land.dim[0]+1),
                    np.linspace(0, mod.land.dim[1], mod.land.dim[1]+1),
                    mod.land[0].rast, cmap='RdBu_r')
ax1 = plt.subplot(gs[0])
mod.plot_phenotype(0, 0, 0)
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('environmental and phenotypic value', rotation=270,
               labelpad=15, y=0.5)
ax2 = plt.subplot(gs[1])
im = plt.pcolormesh(np.linspace(0, mod.land.dim[0], mod.land.dim[0]+1),
                    np.linspace(0, mod.land.dim[1], mod.land.dim[1]+1),
                    mod.land[1].rast, cmap='BrBG_r')
mod.plot_phenotype(0, 1, 1)
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('environmental and phenotypic value', rotation=270,
               labelpad=15, y=0.5)
ax1.set_title('trait_0')
ax2.set_title('trait_1')
plt.show()
fig.tight_layout()
# print out time
print("Model ran in %0.2f seconds." % tot_time)
