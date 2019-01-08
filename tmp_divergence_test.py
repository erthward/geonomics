import geonomics as gnx

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

# set total runtime for each selection coefficient's model
T = 1000
# set selection coefficients
phis = (0.1, 0.05, 0.01)
# set the data directory, and delete it if it already exists (so that we don't
# create mutliple, conflicting sets of data files)
data_dir = './GEONOMICS_mod-divergence_params'
if os.path.isdir(data_dir):
    shutil.rmtree(data_dir)

# create data structures to store allele frequency trajectories and migrations
# for all selection coefficients
allele_freqs = {}
migration_rates = {}

# for each strength of selection, create and run the model, tracking allele
# frequencies in both halves of the environment
for phi in phis:
    # create the model
    mod = gnx.make_model('./test/validation/divergence/divergence_params.py')
    mod.comm[0].gen_arch.traits[0].phi = phi
    # landscape and community will not be randomized between iterations,
    # so I can just extract the non-neutral loci now
    nonneut_loc = mod.comm[0].gen_arch.traits[0].loci

    # burn in the model
    mod.walk(mode='burn', T=200000, verbose=True)

    # create data structures to record allele frequencies in each half of the
    # landscape, and migrations between the halves, at each timestep
    allele_freqs_this_phi = {0: [], 1: []}
    migration_rates_this_phi = {(0, 1): [], (1, 0): []}
    # update it for the starting allele frequencies
    for allele in (0, 1):
        allele_counts = [sum(i.genome[nonneut_loc,
                             :]) for i in mod.comm[0].values(
                             ) if i.e[0] == allele]
        allele_freq = sum(allele_counts) / (2 * len(allele_counts))
        allele_freqs_this_phi[allele].append(allele_freq)
    # walk the model, calculating the 1-allele freq at the non-neutral locus
    # after each timestep, for each half of the landscape
    for t in range(T):
        # get starting environmental values for each individual
        old_env_vals = {i.idx: i.e[0] for i in mod.comm[0].values()}
        mod.walk(1, verbose=True)
        # get allele frequencies for each half of the environment
        for allele in (0, 1):
            allele_counts = [sum(i.genome[nonneut_loc,
                                 :]) for i in mod.comm[0].values(
                                 ) if i.e[0] == allele]
            allele_freq = sum(allele_counts) / (2 * len(allele_counts))
            allele_freqs_this_phi[allele].append(allele_freq)
        # get ending environmental values for each individual
        new_env_vals = {i.idx: i.e[0] for i in mod.comm[0].values()}
        # count all migration events
        mig_0_to_1 = 0
        mig_1_to_0 = 0
        for idx, new_val in new_env_vals.items():
            if idx in [*old_env_vals]:
                if new_val != old_env_vals[idx]:
                    if new_val == 1:
                        mig_0_to_1 += 1
                    else:
                        mig_1_to_0 += 1
        # convert to migration rates
        mig_0_to_1 = mig_0_to_1 / len({i for i in mod.comm[0].values(
                                                            ) if i.e[0] == 1})
        mig_1_to_0 = mig_1_to_0 / len({i for i in mod.comm[0].values(
                                                            ) if i.e[0] == 0})
        # store migration rates
        migration_rates_this_phi[(0, 1)].append(mig_0_to_1)
        migration_rates_this_phi[(1, 0)].append(mig_1_to_0)

    # store allele freq trajectories and migration rates
    # for this selection coefficient
    allele_freqs[phi] = allele_freqs_this_phi
    migration_rates[phi] = migration_rates_this_phi


# define fn to calculate the expected allele-frequency trajectories
# (using migration-selection balance eqxn, pg. 308 Hartl & Clark)
# ('f' refers to allele frequency; 'in' and 'out' refer to freqencies
# inside or outside the population for which delta_q is being calculated)
def mig_sel_bal(s, f0_in, f1_in, f1_out, m_in, m_out, h=0.5):
    delta_q = ((-s * f0_in * f1_in * (f1_in + h * (f0_in - f1_in))) / (
                1 - (s * f1_in * (2 * h * f0_in + f1_in)))) + (
                    m_in * f1_out) - (m_out * f1_in)
    return delta_q


# create data structure to store all expected allele freq trajectories
expected_allele_freqs = {}
# for each phi
for phi in phis:
    # create data structure to store expected allele freq trajectories for
    # each half of the environment for this phi
    expected_this_phi = {0: [], 1: []}
    # for each allele
    for allele in [*expected_this_phi]:
        # generate the 2-tuple keys to use to get in- and out-migration
        # rates for this allele's half of the environment from the
        # migration_rates dict
        m_in_key = (allele, int(allele == 0))
        m_out_key = m_in_key[::-1]
        # append the starting allele freq for this allele's half of the
        # environment as the starting point for the expectation trajectory
        expected_this_phi[allele].append(allele_freqs[phi][allele][0])
        # for each additional timestep
        for t in range(T):
            # get the current 0- and 1-allele freqs for this allele's half
            # of the environment
            curr_f1_in = allele_freqs[phi][allele][t + 1]
            curr_f0_in = 1 - curr_f1_in
            # get the current f1 freq in the other half of the environment
            curr_f1_out = allele_freqs[phi][int(allele == 0)][t + 1]
            # get the previous step's rate of migrations of individuals
            # into this allele's half of the environment
            #curr_m_in = np.mean(migration_rates[phi][m_in_key])
            curr_m_in = migration_rates[phi][m_in_key][t]
            # get the previous step's rate of migration of individuals
            # out of this allele's half of the environment
            #curr_m_out = np.mean(migration_rates[phi][m_out_key])
            curr_m_out = migration_rates[phi][m_out_key][t]
            # TODO: REASON THIS OUT AND FIX COMMENT!
            if allele == 0:
                s = -1*phi #* int(allele == 0)
            else:
                s = phi
            # use those param values to calculated the expected change in
            # freq of the 1-allele in this allele's half of the environment
            expected_delta_q = mig_sel_bal(s=s, f0_in=curr_f0_in,
                                           f1_in=curr_f1_in,
                                           f1_out=curr_f1_out,
                                           m_in=curr_m_in,
                                           m_out=curr_m_out)
            # add the expected change in 1-allele freq to the previous
            # timestep's 1-allele freq to get the expected freq for this
            # timestep
            if allele == 0:
                expected_this_phi[allele].append(
                    expected_this_phi[allele][t] - expected_delta_q)
            else:
                expected_this_phi[allele].append(
                    expected_this_phi[allele][t] - expected_delta_q)
    # store the expected allele frequencies for this phi
    expected_allele_freqs[phi] = expected_this_phi

# plot observed versus expected allele frequencies
fig = plt.figure()
plt.suptitle('Allele-frequency divergence in constrasting environments')
plt.xlabel('time')
plt.ylabel('frequency of 1 allele')
plt.ylim((0, 1))
markers = ['o', 'P', '*']
lines = ['-', '--', ':']
colors = {0: 'blue', 1: 'white'}
line_colors = {0: 'blue', 1: 'black'}
for n, phi in enumerate(phis):
    for allele in (0, 1):
        # plot observed allele frequencies
        plt.plot(range(len(allele_freqs[phi][allele])),
                 allele_freqs[phi][allele], markers[n],
                 color=colors[allele], markeredgecolor='black',
                 markersize=10)
        # plot expected allele frequencies
        plt.plot(range(len(expected_allele_freqs[phi][allele])),
                 expected_allele_freqs[phi][allele], lines[n],
                 color=line_colors[allele], linewidth=1)
plt.show()
