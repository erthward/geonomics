#!/usr/bin/python
# sweep_test.py

# flake8: noqa

import geonomics as gnx

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter as C
from itertools import combinations
import os
from copy import copy
import statsmodels.api as sm
lowess = sm.nonparametric.lowess

# set some plotting params
axlabelsize = 18
titlesize = 20
ticklabelsize=15

#img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper_images/'
#           'final/')
ax_fontdict = {'fontsize': axlabelsize,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': titlesize,
                'name': 'Bitstream Vera Sans'}
color_dict = dict(zip(np.linspace(0,1,3), plt.cm.coolwarm(np.linspace(0,1,3))))


"""
def plot_phenotype_map(ax, mod, title):
    xs = mod.get_x()-0.5
    ys = mod.get_y()-0.5
    zs = mod.get_z(0)
    ax.imshow(mod.land[0].rast, cmap='Wistia')
    if np.all(zs == 1):
        ax.scatter(xs, ys, c=mpl.colors.to_hex(mpl.cm.Wistia_r(256)),
                   edgecolor='black')
    else:
        ax.scatter(xs, ys, c=zs, cmap='Wistia_r', edgecolor='black')
    ax.set_title(title, fontsize=20)
"""


def calc_pi(mod, nonneut_loc, window_width=6):
    print("\nCalculating nucleotide diversity (WILL TAKE A WHILE)...\n")
    # make speciome
    speciome = mod.comm[0]._get_genotypes()
    # create data structure to store pi value for each genomic window
    pi_vals = []
    # for each locus
    for loc in range(speciome.shape[1]):
        # get 5-locus window around the locus
        start = max(loc-window_width, 0)
        # add 6 so that the resulting interval includes 5 loci on either
        # side of focal locus
        stop = min(loc+window_width+1, speciome.shape[1])
        print('\tcalculating for interval [%i, %i)' % (start, stop))
        # get all individs' genotypes for that interval of loci
        interval_idxs = list(range(start, stop))
        # NOTE: get rid of locus 50 itself, if in there, because it is the only
        # non-segregating locus, so it creates an artefactual depression in
        # nucleotide diversity in the area surrounding it
        #if nonneut_loc in interval_idxs:
        #    interval_idxs = [n for n in interval_idxs if n != nonneut_loc]
        intervalome = speciome[:, interval_idxs, :]
        # break into a list by seq (i.e. melt down to 2 seqs per individ)
        intervalome = np.vstack([intervalome[:, :, i] for i in (0, 1)])
        # get freq of each unique sequence
        seq_counts = C([''.join([str(n) for n in seq]) for seq in intervalome])
        seq_freqs = {seq: ct / intervalome.shape[0] for seq,
                     ct in seq_counts.items()}
        # get all pairwise combos of seqs
        combos = [*combinations([*seq_freqs], 2)]
        # get number of nucleotide differences per nucleotide site for each
        # pairwise combo of seqs
        pi_ij = [np.mean([i != j for i, j in zip(c[0], c[1])]) for c in combos]
        # get sum across all pairwise seq combos of mean number of nucleotide
        # differences multiplied by frequencies of each of the two seqs
        pi = 2 * sum([seq_freqs[c[0]] * seq_freqs[c[1]] * pi_ij[n] for n, c in
                     enumerate(combos)])
        pi_vals.append(pi)
    return pi_vals


def calc_mean_fit(mod):
    mean_fit = np.mean(mod.comm[0]._get_fit())
    return mean_fit


def calc_nonneut_loc_freq(mod):
    nonneut_loc_freq = np.mean(np.vstack(
                [i.g[nonneut_loc_idx, :] for i in mod.comm[0].values()]))
    return nonneut_loc_freq


def gen_data_linkage_dist_plot(r2_mat, recomb_rates, L=1001, sweep_loc=500):
    # TODO: Do I need to account for possibility of even number of recombination
    # events between loci when calculating distance? In other words, should
    # I convert the distance to a true cM distance, or just leave as the sum
    # of the interlocus recombination rates?
    r2s = []
    dists = []
    for i in range(L-1):
        for j in range(i, L-1):
            r2 = r2_mat[i, j]
            if not np.isnan(r2):
                dist = np.sum(recomb_rates[i+1:j+1])
                r2s.append(r2)
                dists.append(dist)
    return(r2s, dists)


################################################################

# create data structures to save data to be analyzed/plotted
mean_fit = []
pi = []
ld = []

# build the model
mod = gnx.make_model('./tests/validation/sweep/sweep_params.py')

#get the non-neutral locus
nonneut_loc = mod.comm[0].gen_arch.traits[0].loci[0]
# NOTE: changing this to 0 because only the single non-neutral locus
# will be in each individ's genome array, so idx will be 0
nonneut_loc_idx = 0

# burn the model in
mod.walk(20000, 'burn', True)

# create the figure and its Axes instances
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# calculate LD at outset
print('\ncalculating linkage...\n')
ld.append(gnx.sim.stats._calc_ld(mod.comm[0]))
r2s, dists = gen_data_linkage_dist_plot(ld[0],
                                        mod.comm[0].gen_arch.recombinations._rates)
ax3.plot(dists, r2s, '.r')
ax3.set_xlabel("recombination distance (sum of interlocus recomb. rates)",
              fontdict=ax_fontdict)
ax3.set_ylabel("linkage ($R^{2}$)", fontdict=ax_fontdict)
ax3.tick_params(labelsize=ticklabelsize)

# create a figure to hold the 3 phenotype maps of the sweep
# (one showing the phenotype plot right after fixation,
# the next showing it part-way, and the third after fixation)
map_fig = plt.figure()

# run for 1 main step, so that all individuals have fitness values
mod.walk(1, 'main', True)

# keep reintroducing a mutant and running for 100 timesteps, until
# mutant establishes
sweeping_allele_established = False
lapsed_t = 0
while not sweeping_allele_established:

    mean_fit = [-999]
    pi_at_intro = []

    # store the starting mean fitness
    mean_fit[0] = calc_mean_fit(mod)

    # introduce the adaptive allele in a random individual
    print('\n=================\nINTRODUCING MUTATION\n')
    individ = np.random.choice([*mod.comm[0]])
    chromatid = np.random.binomial(1, 0.5, 1)
    mod.comm[0][individ].g[nonneut_loc_idx, chromatid] = 1

    pi_at_intro.append(calc_pi(mod, nonneut_loc))

    # map the phenotypes
    #map_ax = map_fig.add_subplot(131)
    #map_ax.cla()
    #plot_phenotype_map(map_ax, mod, 'beginning')


    # walk 50 timesteps
    for _ in range(5):
        mod.walk(11)
        # calclate and store mean fitness
        mean_fit.append(calc_mean_fit(mod))

    # check if the allele was established
    nonneut_loc_freq = calc_nonneut_loc_freq(mod)
    sweeping_allele_established = nonneut_loc_freq > 0.05
    if sweeping_allele_established:
        print('\n\nMUTANT ESTABLISHED.\n')
        # calculate and store ending nucleotide diversity
        pi.append(pi_at_intro[-1])


    else:
        print('\n\nMutation lost...\n')

    # increment lapsed_t
    lapsed_t += 50

# plot status after 100 more timesteps
for _ in range(10):
    mod.walk(10)
    # calclate and store mean fitness
    mean_fit.append(calc_mean_fit(mod))
    print('\tallele_freq = %0.3f' % calc_nonneut_loc_freq(mod))

    # map the phenotypes
    #map_ax = map_fig.add_subplot(132)
    #plot_phenotype_map(map_ax, mod, 'during')

# walk until mutant fixed
nonneut_loc_freq = np.mean(np.vstack(
                        [i.g[nonneut_loc_idx, :] for i in mod.comm[0].values()]))
fixed = nonneut_loc_freq == 1
while not fixed:
    mod.walk(10, verbose=True)
    print('\tallele_freq = %0.3f' % calc_nonneut_loc_freq(mod))
    # calclate and store mean fitness
    mean_fit.append(calc_mean_fit(mod))
    # check if fixed
    nonneut_loc_freq = np.mean(np.vstack(
                        [i.g[nonneut_loc_idx, :] for i in mod.comm[0].values()]))
    fixed = nonneut_loc_freq == 1

# calculate and store ending nucleotide diversity
pi.append(calc_pi(mod, nonneut_loc))
print('\ncalculating linkage...\n')

# map phenotype last time
#map_ax = map_fig.add_subplot(133)
#plot_phenotype_map(map_ax, mod, 'after')
#map_fig.show()

# calculate linkage
ld.append(gnx.sim.stats._calc_ld(mod.comm[0]))
r2s, dists = gen_data_linkage_dist_plot(ld[1],
                                        mod.comm[0].gen_arch.recombinations._rates)
ax4.plot(dists, r2s, '.r')
ax4.set_xlabel("recombination distance (sum of interlocus recomb. rates)",
               fontdict=ax_fontdict)
ax4.tick_params(labelsize=ticklabelsize)

# plot stats
ax1.set_title('before sweep', fontdict=ttl_fontdict)
ys = lowess(pi[0], [*range(mod.comm[0].gen_arch.L)], frac=0.01)[:, 1]
ax1.plot(range(mod.comm[0].gen_arch.L), ys, '-k')
ax1.plot([nonneut_loc, nonneut_loc], [0, 1], '--r')
ax1.set_xlabel('genomic_position', fontdict=ax_fontdict)
ax1.set_ylabel('nucleotide diversity', fontdict=ax_fontdict)
ax1.tick_params(labelsize=ticklabelsize)

ax2.set_title('after sweep', fontdict=ttl_fontdict)
ys = lowess(pi[1], [*range(mod.comm[0].gen_arch.L)], frac=0.01)[:, 1]
ax2.plot(range(mod.comm[0].gen_arch.L), ys, '-k')
ax2.plot([nonneut_loc, nonneut_loc], [0, 1], '--r')
ax2.set_xlabel('genomic_position', fontdict=ax_fontdict)
ax2.set_yticks([])

pi_min_lim = 0.95 * min([val for sublist in pi for val in sublist])
pi_max_lim = 1.05 * max([val for sublist in pi for val in sublist])
ax1.set_ylim((pi_min_lim, pi_max_lim))
ax2.set_ylim((pi_min_lim, pi_max_lim))
ax2.tick_params(labelsize=ticklabelsize)

# add vertical space between the first and second rows of plots
fig.subplots_adjust(hspace=.5)

plt.show()
#plt.savefig(os.path.join(img_dir, 'SWEEP_pop_and_nuc_div.pdf'))
print('Sweep complete')

plt.rcParams['figure.figsize'] = [5.5, 4]
fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.plot(np.linspace(0, mod.t, len(mean_fit)), mean_fit, '-', color='#D55E00')
ax.set_xlabel('time', fontdict=ax_fontdict)
ax.set_ylabel('mean fitness', fontdict=ax_fontdict)
ax.set_ylim(((1 - mod.comm[0].gen_arch.traits[0].phi)-0.05, 1.01))
ax.set_xlim((0, mod.t))
ax.set_yticks(np.linspace(0.8, 1.0, 5))
ax.tick_params(labelsize=ticklabelsize)

plt.show()
#plt.savefig(os.path.join(img_dir, 'SWEEP_mean_fit.pdf'))
