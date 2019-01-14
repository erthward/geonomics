#!/usr/bin/python
# yosemite_example.py

import geonomics as gnx
import utils.spatial as spt

import numpy as np
import utils.io as io
import matplotlib.pyplot as plt


# little fn to calculate a habitat raster that's 1 at the center of env-var's
# range and drops toward 0 at extremes
def calc_hab_rast(tmp, l_lim=7, u_lim=12):
    hab = np.ones(tmp.shape)
    hab[tmp < l_lim] = (1 - (tmp[tmp < l_lim] - l_lim) / (tmp.min() - l_lim))
    hab[tmp > u_lim] = (1 - (u_lim - tmp[tmp > u_lim]) / (u_lim - tmp.max()))
    # hab = 1-abs(env_rast-0.5)
    # hab = (hab - hab.min()) / (hab.max() - hab.min())
    return(hab)


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

# burn in, then plot starting population, on both rasters
mod.walk(20000, 'burn')
fig = plt.figure()
ax1 = fig.add_subplot(231)
ax1.set_title('starting population\n(genotypes randomly assigned)')
mod.plot_phenotype(0, 0, 0, size=8)
ax4 = fig.add_subplot(234)
mod.plot_phenotype(0, 0, 1, size=8)

# walk for 500 timesteps, then plot again, before climate change starts
mod.walk(500)
ax2 = fig.add_subplot(232)
ax2.set_title('after 500 timesteps,\n(before climate change begins)')
mod.plot_phenotype(0, 0, 0, size=8)
ax5 = fig.add_subplot(235)
mod.plot_phenotype(0, 0, 1, size=8)

# walk for 1000 more timesteps, then plot again,
# at end of climate-change period
mod.walk(1000)
ax3 = fig.add_subplot(233)
ax3.set_title('after 1500 timesteps\n(at end of period of climate change)')
mod.plot_phenotype(0, 0, 0, size=8)
plt.colorbar()
ax6 = fig.add_subplot(236)
mod.plot_phenotype(0, 0, 1, size=8)
plt.colorbar()
fig.suptitle(('Evolutionary response of 100-locus additive trait to climate '
              'change in Yosemite region, N = ~14300 individuals\n'
              'row 1: temperature rasters (i.e. selective environment); '
              'row 2: habitat quality rasters (i.e. carrying capacity)'))
ax1.set_ylabel('temperature rasters')
ax4.set_ylabel('habitat rasters')

# TODO: 
    # add colorbars for phenotype and for rasters
    # create some interpolated/rasterized map of mean phenotye 



plt.show()

# TODO run analyses?
