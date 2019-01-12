#!/usr/bin/python
# yosemite_example.py

import geonomics as gnx
import utils.io as io
import matplotlib.pyplot as plt


# little fn to calculate a habitat raster that's 1 at the center of env-var's
# range and drops toward 0 at extremes
def calc_hab_rast(env_rast):
    hab = 1-abs(env_rast-0.5)
    hab = (hab - hab.min()) / (hab.max() - hab.min())
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
# fut_tmp, fut_min, fut_max = spt.scale_raster(fut_tmp, -1.37, 19.11)

# add the future tmp to params
params.landscape.layers['tmp'].change[0]['change_rast'] = fut_tmp

# set the scale-min and scale-max values for the tmp layer
scale_min = min(tmp.min(), fut_tmp.min())
scale_max = min(tmp.max(), fut_tmp.max())
params.landscape.layers['tmp'].init.file['scale_min_val'] = scale_min
params.landscape.layers['tmp'].init.file['scale_max_val'] = scale_max

# use tmp rasters to create start and end habitat rasters, and add to params
hab = calc_hab_rast(tmp)
fut_hab = calc_hab_rast(fut_tmp)
# and add them to the params
params.landscape.layers['hab'].init.defined['rast'] = hab
params.landscape.layers['hab'].change[0]['change_rast'] = fut_hab

# create the model
mod = gnx.make_model(params)

# run the model
print('RUNNING MODEL...\n\n')
# mod.run(verbose=True)
mod.walk(20000, 'burn')
fig = plt.figure()
mod.plot_phenotype(0, 0, 0)
mod.walk(2000)
fig = plt.figure()
mod.plot_phenotype(0, 0, 0)

# run analyses for all 3 iterations
