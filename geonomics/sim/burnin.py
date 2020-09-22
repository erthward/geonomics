#!/usr/bin/python
#burn_in.py

'''
Functions for testing burn-in stationarity.
'''

# other imports
import numpy as np
from numpy import random as r
from matplotlib import colors
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import cKDTree
from collections import Counter as C
from scipy.stats import ttest_rel as ttest
from statsmodels.tsa.stattools import adfuller as adf
from sklearn.preprocessing import normalize


class SpatialTester:
    def __init__(self, spp):
        # get the land dimensions
        self.dim = spp._land_dim

        # create the initial counts and diff arrays
        self.counts = np.zeros(self.dim)
        self.diff = np.zeros(self.dim)

        # create the stats dict
        self.stats = {np.mean: [],
                      #np.median: [],
                      np.std: []
                     }

        # update the counts and stats
        self.update(spp)


    def update_counts_cell(self, i, j, v):
        self.counts[i,j] = v


    def update(self, spp):
        # get a copy of the previous counts
        copy = np.copy(self.counts)
        # get the species, then get a counter of number of individuals
        # within each cell
        counter = C([(int(x), int(y)) for x,y in zip(spp._get_x(),
                                                     spp._get_y())])
        # convert that counter into an array of cell counts
        for i in range(self.dim[0]):
            for j in range(self.dim[1]):
                self.update_counts_cell(i, j, counter.get((j,i), 0))
        # update diff
        self.diff = self.counts - copy
        # update stats
        for fn, v in self.stats.items():
            self.stats[fn].append(fn(self.diff))


    def plot_diffs(self, spp, bins=50, alpha=0.5):
        fig = plt.figure()
        # histogram
        ax1 = fig.add_subplot(121)
        fig.suptitle('Histogram of individual-count changes by landscape cell')
        plt.hist(self.diff, bins=bins, alpha=alpha)
        # map
        ax2 = fig.add_subplot(122)
        fig.suptitle('Map of individual-count changes')
        colors.DivergingNorm(vcenter=0.)
        plt.imshow(self.diff, origin='lower', cmap='coolwarm')
        plt.colorbar()
        plt.scatter(spp._get_x() - 0.5, spp._get_y() - 0.5, c='black', s=15)
        plt.show()

    def run_test(self, num_timesteps_back, alpha=0.05):
        results = []
        for fn, data in self.stats.items():
            try:
                adf_res = adf(data[-num_timesteps_back:])[1] < alpha
            except ValueError as e:
                adf_res = None
            try:
                ttest_res = ttest(data[int(-num_timesteps_back):
                                       int(-num_timesteps_back/2)],
                              data[int(-num_timesteps_back/2):])[1] > alpha
            except ValueError as e:
                ttest_res = None
            results.append(adf_res and ttest_res)
        return np.all(results)


def _test_adf_threshold(spp, num_timesteps_back, alpha=0.05):
    result = adf(spp.Nt[-num_timesteps_back:])[1] < alpha
    return result


def _test_t_threshold(spp, num_timesteps_back, alpha=0.05):
    num_timesteps_back += num_timesteps_back % 2
    result = ttest(spp.Nt[int(-num_timesteps_back): int(-num_timesteps_back/2)],
                                spp.Nt[int(-num_timesteps_back/2):])[1] > alpha
    return result
