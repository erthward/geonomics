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


# number its to run for each recomb rate scenario and program (msprime, gnx)
n_its = 2

# effective population size for the sims
Ne = 10

# read the parameters files for the fixed and custom recombination models
params_filename = ('/home/deth/Desktop/stuff/berk/research/'
                   'projects/sim/geonomics/tests/validation/'
                   'recomb/validate_recomb_params_recomb_%s.py')
params_fix = gnx.read_parameters_file(params_filename % 'fixed')
params_cus = gnx.read_parameters_file(params_filename % 'custom')


##################################
# run the fixed-recombination sims

# recombination rates and sequence lengths for fixed-rate simulations
rates_and_lengths = [[1e-2, 1e3]]

# struct to hold the numbers of trees from both programs for all recomb rates
num_trees_dict_fix = {'msp': {},
                      'gnx': {}
                 }
"""
for rate, length in rates_and_lengths:

    print('\nRUNNING SIMS FOR RATE %.2E, LENGTH %i\n' % (rate, int(length)))

    # structs to hold each it's num trees for each program
    msp_n_trees = []
    gnx_n_trees = []

    # make the fixed-recombination params for this rate-length combo
    params = deepcopy(params_fix)
    params['comm']['species']['spp_0']['gen_arch']['L'] = int(length)
    params['comm']['species']['spp_0']['gen_arch']['r_distr_alpha'] = rate

    for i in range(n_its):
        print('\n\tITERATION: %i\n' % i)
        # run the msprime simulation and get its number of trees
        print('\n\t\tMSPRIME RUNNING...\n')
        ts = msprime.simulate(100,
                              Ne=Ne,
                              length=length,
                              recombination_rate=rate)
        msp_n_trees.append(ts.num_trees)

        # get the max root time across the tree sequence,
        # then round up and use this to run a gnx simulation
        max_t = int(np.ceil(ts.max_root_time))

        # run a geonomics model for that length of time
        print('\n\t\tGNX RUNNING...\n')
        #assert True == False
        mod = gnx.make_model(params)
        mod.walk(100000, 'burn')
        for t in range(max_t):
            mod.walk(1)

        # get the gnx model's number of trees
        mod.comm[0]._sort_simplify_table_collection()
        ts = mod.comm[0]._tc.tree_sequence()
        gnx_n_trees.append(ts.num_trees)

    # add this rate's n_trees lists to the overall data structure
    num_trees_dict_fix['msp'][rate] = msp_n_trees
    num_trees_dict_fix['gnx'][rate] = gnx_n_trees

# plot hists
fix_fig = plt.figure()
ax_msp = fix_fig.add_subplot(121)
ax_gnx = fix_fig.add_subplot(122)
ax_msp.set_title('msprime')
ax_gnx.set_title('geonomics')
axes = [ax_msp, ax_gnx]
for ax in axes:
    ax.set_xlabel('number of trees')
    ax.set_ylabel('frequency')
for rate, length in rates_and_lengths:
    ax_msp.hist(num_trees_dict_fix['msp'][rate],
                bins=25, alpha=0.5,
                label='rate: %.2E; length: %.2E' % (rate, length))
    ax_gnx.hist(num_trees_dict_fix['gnx'][rate],
                bins=25, alpha=0.5,
                label='rate: %.2E; length: %.2E' % (rate, length))
ax_msp.legend(title='recomb. rate')
ax_gnx.legend(title='recomb. rate')
fix_fig.show()
"""

###################################
# run the custom-recombination sims

# read the custom gen_arch file and create an msprime recomb-map object from it
gen_arch_file = (params_filename % 'custom').rstrip(
                    '.py') + '_spp-0_gen_arch.csv'
df = pd.read_csv(gen_arch_file)
loci = [*df.locus.values]
rates = [*df.r.values[1:]] + [0]
recomb_map = msprime.RecombinationMap(loci, rates, num_loci=len(loci))

# run the msprime simulation and get its number of trees
print('\n\t\tMSPRIME RUNNING...\n')
ts = msprime.simulate(100,
                      Ne=1000,
                      recombination_map=recomb_map)

# code taken and tweaked from Jerome Kelleher and Konrad Lohse's
# chapter on msprime in Statistical Population Genomics textbook (2020);
# GitHub repo: https://github.com/StatisticalPopulationGenomics/msprime
msp_breakpoints = np.array(list(ts.breakpoints()))
# Now we get the positions and rates from the recombination
# map and plot these using 100 bins
positions = np.array(recomb_map.get_positions()[1:])-0.5
rates = np.array(recomb_map.get_rates()[:-1])
num_bins = 100
v, bin_edges, _ = scipy.stats.binned_statistic(positions, rates,
                                               bins=num_bins)
x = bin_edges[:-1][np.logical_not(np.isnan(v))]
y = v[np.logical_not(np.isnan(v))]
fig_cus, ax_true = plt.subplots()
ax_true.plot(x, y, ':k', label='true rate')
ax_true.set_ylabel("recombination rate", fontsize=18)
ax_true.set_xlabel("locus", fontsize=18)
ax_true.tick_params(labelsize=15)
# then plot the observed breakpoint densities from msprime
ax_obs = ax_true.twinx()
v, bin_edges = np.histogram(msp_breakpoints, num_bins, density=True)
ax_obs.plot(bin_edges[:-1], v, '-b', label='msprime')
ax_obs.set_ylabel("breakpoint density", fontsize=18)
ax_obs.tick_params(labelsize=15)

#num_trees_dict_cus['msp'].append(ts.num_trees)

# run a geonomics model for that length of time
print('\n\t\tGNX RUNNING...\n')
#assert True == False
mod = gnx.make_model(params_cus)
mod.walk(100000, 'burn')
mod.walk(20)

# get the gnx model's number of trees
ts_gnx = mod.get_tree_sequence(0)

# add the gnx breakpoint densities to the plot
gnx_breakpoints = np.array(list(ts_gnx.breakpoints()))
v, bin_edges = np.histogram(gnx_breakpoints, num_bins, density=True)
ax_obs.plot(bin_edges[:-1], v, '-r', label='geonomics')

ax_obs.legend()

fig_cus.show()
