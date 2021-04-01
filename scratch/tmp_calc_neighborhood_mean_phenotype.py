#!/usr/bin/python

#####################
#TODO
#1.  Add a "timing" argument, that I can set to false to run the model without
#plotting, to get an accurate assessment of run time

#2. Debug all the plotting stuff
#####################

# geonomics imports
from geonomics.utils.viz import _check_display

# other imports
import os
import numpy as np
import matplotlib as mpl
_check_display()
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d.axes3d import Axes3D
from nlmpy import nlmpy
import time
import rasterio as rio
from scipy.stats import norm


gnx_dir = os.path.split(__file__)[0]
gnx_dir = os.path.join(*os.path.split(gnx_dir)[:-1])
DATA_PATH = os.path.join(gnx_dir, "data", "yosemite_demo")


def calc_neighborhood_mean_phenotype(dim, comm, window_width=8,
                                     rel_bandwidth=0.5):
    bandwidth = window_width * rel_bandwidth
    print(bandwidth)
    # array to store each cell's mean phenotype value
    mean_z = np.ones(dim)
    # calc half-window_width
    hww = int(window_width / 2)
    # create the kernel to use for weighting individuals' phenotypes
    pd = norm(0, bandwidth)
    # calculate the max prob dens (to normalize weights in the weights lists
    # below)
    max_weight = pd.pdf(0)
    # loop over cells
    for i in range(dim[1]):
        for j in range(dim[0]):
            # get window around each cell
            i_min = max(i + 0.5 - hww, 0)
            i_max = min(i + 0.5 + hww, dim[1])
            j_min = max(j + 0.5 - hww, 0)
            j_max = min(j + 0.5 + hww, dim[0])
            # get all phenotypes in the window, Gaussian-weighted by their
            # distances from the window center
            zs_in_window = []
            dists_from_center = []
            effective_zs = []
            for ind in comm.values():
                if ((i_min <= ind.y <= i_max)
                    and (j_min <= ind.x <= j_max)):
                    zs_in_window.append(ind.z)
                    dist = np.sqrt((i + 0.5 - ind.x)**2 + (j + 0.5 - ind.y)**2)
                    dists_from_center.append(dist)
            # average the window's phenotypes and add to mean_z
            # NOTE: if there are no individuals in the window then a NaN is
            # returned
            weights = [pd.pdf(d)/max_weight for d in dists_from_center]
            print(weights)
            effective_zs = [zs_in_window[n] * weights[n] for n in range(
                                                                len(weights))]
            mean_eff_zs = np.mean(effective_zs)
            if np.isnan(mean_eff_zs):
                mean_eff_zs = 0
            mean_z[i, j] = np.mean(effective_zs)

    return mean_z

