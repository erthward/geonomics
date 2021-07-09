#!/usr/bin/python

# flake8: noqa

# tmp_divergence_test.py

import geonomics as gnx

import os
import numpy as np
import shutil
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.api as sm

# set some plotting params
#img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
#           'img/final/')
titlesize=20
axlabelsize=18
ticklabelsize=15
ax_fontdict = {'fontsize': axlabelsize,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': titlesize,
                'name': 'Bitstream Vera Sans'}


# define a function to calculate the mean difference between phenotype
# and environment for a species
def calc_mean_z_e_diff(spp, trait_num=0):
    zs = spp._get_z().ravel()
    es = spp._get_e(lyr_num=spp.gen_arch.traits[trait_num].lyr_num)
    mean_diff = np.mean(np.abs(zs - es))
    return mean_diff


# set total runtime for each selection coefficient's model
T = 1000
# set selection coefficients
phis = (0.1, 0.05, 0.01)
# set the data directory, and delete it if it already exists (so that we don't
# create mutliple, conflicting sets of data files)
data_dir = './GNX_mod-divergence_params'
if os.path.isdir(data_dir):
    shutil.rmtree(data_dir)

# create data structures to store allele frequency trajectories,migrations,
# and phenotype-environment diff values for all selection coefficients
allele_freqs = {}
migration_rates = {}
z_e_diffs = {}

# create a dict to store the start and end phenotype-environment (z-e)
# correlations for each phi
z_e_corrs = {}

# for each strength of selection, create and run the model, tracking allele
# frequencies in both halves of the environment
for phi in phis:
    # create the model
    mod = gnx.make_model(('./tests/validation/divergence/'
                         'divergence_params.py'))
    mod.comm[0].gen_arch.traits[0].phi = phi
    # landscape and community will not be randomized between iterations,
    # so I can just extract the non-neutral loci now
    nonneut_loc = mod.comm[0].gen_arch.traits[0].loci

    # burn in the model
    mod.walk(mode='burn', T=200000, verbose=True)

    # using tskit?
    use_tskit = mod.comm[0].gen_arch.use_tskit

    # create data structures to record allele frequencies in each half of the
    # landscape, migrations between the halves, and z-e diffs at each timestep
    allele_freqs_this_phi = {0: [], 1: []}
    migration_rates_this_phi = {(0, 1): [], (1, 0): []}
    z_e_diffs_this_phi = {0: calc_mean_z_e_diff(mod.comm[0])}

    # run and store the starting z-e correlation
    z_e_corr_dict = {}
    z_e_corr_dict[0] = sm.Logit(mod.comm[0]._get_z(),
                                mod.comm[0]._get_e(0).reshape(
                                                (len(mod.comm[0]), 1))).fit()

    # update it for the starting allele frequencies
    for allele in (0, 1):
        if use_tskit:
            allele_counts = [np.sum(i.g == allele) for i in mod.comm[0].values(
                             ) if i.e[0] == allele]
        else:
            nonneut_loc = mod.comm[0].gen_arch.traits[0].loci[0]
            allele_counts = [np.sum(i.g[nonneut_loc,
                :] == allele) for i in mod.comm[0].values() if i.e[0] == allele]
        allele_freq = sum(allele_counts) / (2 * len(allele_counts))
        # flip the freq, if tracking the 0 allele
        if allele == 0:
            allele_freq = 1 - allele_freq
        allele_freqs_this_phi[allele].append(allele_freq)
    # walk the model, calculating the 1-allele freq at the non-neutral locus
    # after each timestep, for each half of the landscape
    for t in range(T):
        # get starting environmental values for each individual
        old_env_vals = {i.idx: i.e[0] for i in mod.comm[0].values()}
        mod.walk(1, verbose=True)
        # get allele frequencies for each half of the environment
        for allele in (0, 1):
            if use_tskit:
                allele_counts = [np.sum(i.g == allele) for i in mod.comm[
                    0].values() if i.e[0] == allele]
            else:
                nonneut_loc = mod.comm[0].gen_arch.traits[0].loci[0]
                allele_counts = [np.sum(i.g[nonneut_loc,
                                :] == allele) for i in mod.comm[
                                0].values() if i.e[0] == allele]
            allele_freq = sum(allele_counts) / (2 * len(allele_counts))
            # flip the freq, if tracking the 0 allele
            if allele == 0:
                allele_freq = 1 - allele_freq
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

        # store new z-e diffs
        z_e_diffs_this_phi[t] = calc_mean_z_e_diff(mod.comm[0])

    # store allele freq trajectories and migration rates
    # for this selection coefficient
    allele_freqs[phi] = allele_freqs_this_phi
    migration_rates[phi] = migration_rates_this_phi
    z_e_diffs[phi] = z_e_diffs_this_phi
    z_e_corr_dict[1] = sm.Logit(mod.comm[0]._get_z(),
                                mod.comm[0]._get_e(0).reshape(
                                                (len(mod.comm[0]), 1))).fit()
    z_e_corrs[phi] = z_e_corr_dict


# define fn to calculate the expected allele-frequency trajectories
# (using migration-selection balance eqxn, Eqxn: 6.25, pg. 308 Hartl & Clark)
# (p and q are 0- and 1-allele frequencies in the focal population,
# i.e. the population for which delta_q is being calculated;
# q_star is the 1-allele frequency in the other population)
# NOTE: assumes HWE, rarity of recessive allele, and 'sufficiently weak'
# selection
def mig_sel_bal(s, p, q, q_star, m_in, m_out, h=0.5):
    delta_q = ((-1*s * p * q * (q + (h * (p - q)))) / (
                1 - (s * q * (2 * h * p + q)))) + (
                    m_in * q_star) - (m_out * q)
    return delta_q


# create data structure to store all expected allele freq trajectories
expected_allele_freqs = {}
# for each phi
for phi in phis:
    # create data structure to store expected allele freq trajectories for
    # each half of the environment for this phi
    expected_this_phi = {0: [], 1: []}
    # append the starting allele freq for each allele's half of the
    # environment as the starting point for the expectation trajectory
    [expected_this_phi[allele].append(
            allele_freqs[phi][allele][0]) for allele in (0, 1)]
    for t in range(T):
        for allele in [*expected_this_phi]:
            # get the current 0- and 1-allele freqs for this allele's half
            # of the environment
            curr_f1 = expected_this_phi[allele][t]
            curr_f0 = 1 - curr_f1
            # get the current f1 freq in the other half of the environment
            curr_f1_star = expected_this_phi[int(allele == 0)][t]
            curr_f0_star = 1 - curr_f1_star
            # get the previous step's rate of migrations of individuals
            # into this allele's half of the environment
            curr_m_1_0 = migration_rates[phi][(1, 0)][t]
            # get the previous step's rate of migration of individuals
            # out of this allele's half of the environment
            curr_m_0_1 = migration_rates[phi][(0, 1)][t]
            params = [curr_f1, curr_f0, curr_f1_star, curr_m_1_0,
                      curr_m_0_1, phi]
            for param in params:
                assert 0 <= param <= 1, ("Param outside the 0-1 range! Here"
                                         " are all params: %s") % str(params)
            # set selection coefficient and calculate expected allele freq
            # changes and resulting 1-allele frequencies correctly
            if allele == 0:
                expected_delta_q = mig_sel_bal(s=phi, p=curr_f0,
                                               q=curr_f1,
                                               q_star=curr_f1_star,
                                               m_in=curr_m_1_0,
                                               m_out=curr_m_0_1)
                expected_delta_1 = expected_delta_q
            else:
                expected_delta_q = mig_sel_bal(s=phi, p=curr_f1,
                                               q=curr_f0,
                                               q_star=curr_f0_star,
                                               m_in=curr_m_0_1,
                                               m_out=curr_m_1_0)
                expected_delta_1 = -1*expected_delta_q
            # add the expected change in 1-allele freq to the previous
            # timestep's 1-allele freq to get the expected freq for this
            # timestep
            next_expected = expected_this_phi[allele][t] + expected_delta_1
            next_expected = np.clip(next_expected, 0, 1)
            expected_this_phi[allele].append(next_expected)
    # store the expected allele frequencies for this phi
    expected_allele_freqs[phi] = expected_this_phi

# plot observed versus expected allele frequencies for all phi values
plt.rcParams['figure.figsize'] = [8, 6]
fig = plt.figure()
ax = fig.add_subplot(111)
# ax.set_title(('Allele-frequency divergence in constrasting environments;'
#              '\nobserved (markers) versus predicted (lines);'
#              '%i timesteps, ~1400 individuals') % T)
plt.xlabel('time', fontdict=ax_fontdict)
plt.ylabel('frequency of \'1\' allele', fontdict=ax_fontdict)
plt.ylim((0, 1))
ax.tick_params(labelsize=ticklabelsize)
#markers = ['o', 'P', '*']
#lines = ['-', '--', ':']
cmap = plt.cm.RdBu_r
#colors = {0: cmap(40),
#          1: cmap(255 - 40)}
#line_colors = {**colors}
colors = {0: {0.01: '#c0d6fc',
              0.05: '#5c8ee6',
              0.1:  '#072d70'},
          1: {0.01: '#ffa8a8',
              0.05: '#eb5757',
              0.1:  '#940606'}
         }
plotted_line_list = []

for n, phi in enumerate(phis):
    for allele in (0, 1):
        # plot observed allele frequencies
        ax.plot(range(len(allele_freqs[phi][allele])),
                allele_freqs[phi][allele], '-',
                color=colors[allele][phi],
                linewidth=2) 
        line, = plt.plot(range(len(expected_allele_freqs[phi][allele])),
                         expected_allele_freqs[phi][allele], ':',
                         color=colors[allele][phi], linewidth=2)
        plotted_line_list.append(line)

#labels=['allele %i: phi = %0.2f' %  (n, phi % 2) for phi in phis for n in (0, 1)]
#labels = [labels[int(i/2)] if i%2 == 0 else '' for i in range(len(labels)*2)]
#ax.legend(labels=labels, loc=(1.02, 0.3), fontsize=15)

phi_colors = ['#3d3d3d', '#808080', '#c4c4c4']
custom_lines = [Line2D([0], [0], color=colors[allele][0.05],
                        lw=2) for allele in (0, 1)]
custom_labels = ['$%i$' % allele for allele in (0, 1)]
leg1 = ax.legend(custom_lines, custom_labels, loc=(1.02, 0.8), fontsize=15,
          title='Allele')
leg1.get_title().set_fontsize('17')
custom_lines = [Line2D([0], [0], color=col,
                       lw=2) for col in phi_colors]
custom_labels = ['$\phi=%0.2f$' % phi for phi in phis]
leg2 = ax.legend(custom_lines, custom_labels, loc=(1.02, 0.4), fontsize=15,
          title='Strength of\nSelection')
leg2.get_title().set_fontsize('17')
custom_lines = [Line2D([0], [0], color='black',
                       linestyle=style) for style in ('-', ':')]
custom_labels = ['observed', 'expected']
leg3 = ax.legend(custom_lines, custom_labels, loc=(1.02, 0.15), fontsize=15,
          title='Data')
leg3.get_title().set_fontsize('17')
ax.add_artist(leg1)
ax.add_artist(leg2)



#plt.savefig(os.path.join(img_dir, 'DIVERGENCE_allele_freqs.pdf'))

plt.rcParams['figure.figsize'] = [6, 6]
fig2 = plt.figure()
# ax2.set_title(('Population, colored by phenotype (outer circle) '
#              'and fitness (inner circle)'))
mod.plot_fitness(0, 0, 0, fitness_cbar=False)
plt.show()
#plt.savefig(os.path.join(img_dir, 'DIVERGENCE_pop_plot.pdf'))

# plot z-e diffs
fig3 = plt.figure()
z_e_ax = fig3.add_subplot(111)

for n, phi in enumerate(phis):
    diffs = z_e_diffs[phi]
    plt.plot(range(len(diffs)), [*diffs.values()], color=phi_colors[n])

z_e_ax.set_xlabel('time', fontdict=ax_fontdict)
z_e_ax.set_ylabel('mean phenotype-environment difference', fontdict=ax_fontdict)
z_e_ax.legend(labels=['$\phi$ = %0.3f' % phi for phi in phis[::-1]],
              loc='best', fontsize=15)
z_e_ax.tick_params(labelsize=ticklabelsize)
plt.subplots_adjust(left=0.05, bottom=0.085, right=0.985, top=0.985,
                       wspace=0, hspace=0)
plt.show()


# print out before and after z-e regression results
print(('\n\nResults of logistic regressions of individual phenotypes on '
       'individuals\' environmental values:\n'))
for phi in phis[::-1]:
    print('\tPHI = %0.2f:' % phi)
    print('\t----------')
    print('\tpseu.R^2\tp-val')
    print('start:  %0.3f\t\t%0.3f' % (z_e_corrs[phi][0].prsquared,
                                      z_e_corrs[phi][0].pvalues[0]))
    print('end:    %0.3f\t\t%0.3f' % (z_e_corrs[phi][1].prsquared,
                                      z_e_corrs[phi][1].pvalues[0]))
    print('\n\n')
