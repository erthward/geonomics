"""
Validate Geonomics' recombination model against msprime's, as suggested by
review #2 of the MBE methods paper
"""

import geonomics as gnx
import msprime

import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy


# effective population size for the msprime sim
Ne = 1000

# number of bins to use to summarize recombination rates
num_bins = 250

# read the parameters file
params_filename = ('/home/deth/Desktop/UCB/research/'
                   'projects/sim/geonomics/tests/validation/'
                   'recomb/recomb_params.py')
params = gnx.read_parameters_file(params_filename)

# read the custom gen_arch file and create an msprime recomb-map object from it
gen_arch_file = params_filename.rstrip('.py') + '_spp-0_gen_arch.csv'
df = pd.read_csv(gen_arch_file)
loci = [*df.locus.values]
rates = [*df.r.values[1:]] + [0]
recomb_map = msprime.RecombinationMap(loci, rates, num_loci=len(loci))

# run the msprime simulation and get its number of trees
print('\n\t\tMSPRIME RUNNING...\n')
ts = msprime.simulate(1000,
                      Ne=Ne,
                      recombination_map=recomb_map)

# code taken and tweaked from Jerome Kelleher and Konrad Lohse's
# chapter on msprime in Statistical Population Genomics textbook (2020);
# GitHub repo: https://github.com/StatisticalPopulationGenomics/msprime
msp_breakpoints = np.array(list(ts.breakpoints()))
# Now we get the positions and rates from the recombination
# map and plot these using 100 bins
positions = np.array(recomb_map.get_positions()[1:])-0.5
rates = np.array(recomb_map.get_rates()[:-1])
num_bins = num_bins
v, bin_edges, _ = scipy.stats.binned_statistic(positions, rates,
                                               bins=num_bins)
x = bin_edges[:-1][np.logical_not(np.isnan(v))]
y = v[np.logical_not(np.isnan(v))]
fig_cus, ax_true = plt.subplots()
ax_true.plot(x, y, ':k')#, label='true rate')
ax_true.set_ylabel("mean recombination rate", fontsize=18)
ax_true.set_xlabel("locus", fontsize=18)
ax_true.tick_params(labelsize=15)
# then plot the observed breakpoint densities from msprime
ax_obs = ax_true.twinx()
v, bin_edges = np.histogram(msp_breakpoints, num_bins, density=True)
ax_obs.plot(bin_edges[:-1], v, '-b', label='msprime')
ax_obs.set_ylabel("breakpoint density", fontsize=18)
ax_obs.tick_params(labelsize=15)

# run a geonomics model for that length of time
print('\n\t\tGNX RUNNING...\n')
#assert True == False
mod = gnx.make_model(params)
mod.walk(100000, 'burn')
mod.walk(20)

# get the gnx model's number of trees
ts_gnx = mod.get_tree_sequence()

# add the gnx breakpoint densities to the plot
gnx_breakpoints = np.array(list(ts_gnx.breakpoints()))
v_gnx, bin_edges_gnx = np.histogram(gnx_breakpoints, num_bins, density=True)
ax_obs.plot(bin_edges_gnx[:-1], v_gnx, '-r', label='geonomics: tskit')
ax_obs.legend(prop={"size":16})


# run the same model, without using tskit
params['comm']['species']['spp_0']['gen_arch']['use_tskit'] = False
mod = gnx.make_model(params)
mod.walk(100000, 'burn')
mod.walk(20)

# function to recover breakpoints from subsetters
def get_subsetter_bps(subsetter):
    path = np.array([int(n) for n in sub][::2])
    bps = np.where((path[:-1] - path[1:]) != 0)[0] + 0.5
    #jitter = np.random.normal(0, 0.0000001, len(bps))
    #bps = bps + jitter
    return bps

gnx_notskit_breakpoints = []
for sub in mod.comm[0].gen_arch.recombinations._subsetters.values():
    gnx_notskit_breakpoints.extend(get_subsetter_bps(sub))
v_gnx_nots, bin_edges_gnx_nots = np.histogram(gnx_notskit_breakpoints,
                                    num_bins, density=True)
ax_obs.plot(bin_edges_gnx_nots[:-1], v_gnx_nots,
            '-y', label='geonomics: no tskit')
ax_obs.legend(prop={"size":16})


fig_cus.show()
