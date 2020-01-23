#!/usr/bin/python
# sweep_test.py

# flake8: noqa

import geonomics as gnx

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter as C
from itertools import combinations
import os
from copy import copy

# set some plotting params
img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
           'img/final/')
ax_fontdict = {'fontsize': 12,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': 14,
                'name': 'Bitstream Vera Sans'}
color_dict = dict(zip(np.linspace(0,1,3), plt.cm.coolwarm(np.linspace(0,1,3))))

def calc_pi(mod, nonneut_loc, window_width=6):
    print("\nCalculating nucleotide diversity (WILL TAKE A WHILE)...\n")
    # make speciome
    speciome = np.stack([i.g for i in mod.comm[0].values()])
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
                [i.g[nonneut_loc, :] for i in mod.comm[0].values()]))
    return nonneut_loc_freq


# function for plotting species colored by phenotype, since I don't yet know
# why the landscape color seems to be incorrectly inverted (i.e. landscape
# of all 1s is plotting as all blue, not all red)
def plot_by_phenotype(mod):
    cols = [color_dict[np.mean(ind.g[mod.comm[0].gen_arch.traits[0].loci[0],
                                        :])] for ind in mod.comm[0].values()]
    y_cell_bds, x_cell_bds = [np.linspace(0, dim,
                                    dim+1) for dim in mod.land[0].rast.shape]
    plt.pcolormesh(x_cell_bds, y_cell_bds, mod.land[0].rast, cmap='coolwarm')
    #plt.imshow(mod.land[0].rast, cmap='RdBu')
    plt.scatter(mod.comm[0]._get_x(),# - 0.5,
                mod.comm[0]._get_y(),# - 0.5,
                c=cols, cmap='coolwarm', edgecolors='black', s=10)
    plt.axis('scaled')
    plt.xticks([])
    plt.yticks([])


# create data structures to save data to be analyzed/plotted
mean_fit = []
pi = []
ld = []

# build the model
mod = gnx.make_model('./tests/validation/sweep/sweep_params.py')
# set the non-neutral locus to the central locus (i.e. 50, in a 101-length
# genomic arch.), then set the starting 1-allele freq of the
# non-neutral allele to 0
#nonneut_loc = 50
#mod.comm[0].gen_arch.traits[0].loci[0] = nonneut_loc
#mod.comm[0].gen_arch.p[nonneut_loc] = 0

#get the non-neutral locus
nonneut_loc = mod.comm[0].gen_arch.traits[0].loci[0]

# burn the model in, then run for 50 timesteps (to make sure pop size is
# at equilibrium)
mod.walk(20000, 'burn', True)
mod.walk(50, 'main', True)

# create the figure and first Axes instance
fig = plt.figure()
ax1 = fig.add_subplot(241)
ax1.set_title('Population at beginning of sweep\n(colored by phenotype)')

# calculate LD at outset
ld.append(gnx.sim.stats._calc_ld(mod.comm[0]))

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
    mod.comm[0][individ].g[nonneut_loc, chromatid] = 1

    pi_at_intro.append(calc_pi(mod, nonneut_loc))

    # plot the beginning of the sweep, with the mutant individual
    plt.cla()
    plot_by_phenotype(mod)
    ax1.set_title('just after mutation', fontdict=ttl_fontdict)

    # walk 50 timesteps
    for _ in range(5):
        mod.walk(10)
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
# calculate and store nucleotide diversity again
pi.append(calc_pi(mod, nonneut_loc))
ax2 = fig.add_subplot(242)
plt.title('200 timesteps into sweep',  fontdict=ttl_fontdict)
plot_by_phenotype(mod)

# walk until mutant fixed
nonneut_loc_freq = np.mean(np.vstack(
                        [i.g[nonneut_loc, :] for i in mod.comm[0].values()]))
fixed = nonneut_loc_freq == 1
while not fixed:
    mod.walk(10, verbose=True)
    print('\tallele_freq = %0.3f' % calc_nonneut_loc_freq(mod))
    # calclate and store mean fitness
    mean_fit.append(calc_mean_fit(mod))
    # check if fixed
    nonneut_loc_freq = np.mean(np.vstack(
                        [i.g[nonneut_loc, :] for i in mod.comm[0].values()]))
    fixed = nonneut_loc_freq == 1

# plot status after sweep is complete
ax3 = fig.add_subplot(243)
ax3.set_title(('just after sweep (t = %i)') % (mod.t - lapsed_t + 50 + 50),
              fontdict=ttl_fontdict)
plot_by_phenotype(mod)
# calculate and store ending nucleotide diversity
pi.append(calc_pi(mod, nonneut_loc))

# walk model 5000 more timesteps, then plot status
mod.walk(5000, verbose=True)
ax4 = fig.add_subplot(244)
ax4.set_title('5000 timesteps after sweep (t = %i)' % (mod.t - lapsed_t + 50 + 50),
              fontdict=ttl_fontdict)
plot_by_phenotype(mod)

# calculate and store nucleotide diversity again
pi.append(calc_pi(mod, nonneut_loc))
# and calculate and store LD again
ld.append(gnx.sim.stats._calc_ld(mod.comm[0]))

# plot stats
ax5 = fig.add_subplot(245)
ax5.plot(range(mod.comm[0].gen_arch.L), pi[0])
ax5.plot([nonneut_loc, nonneut_loc], [0, 1], '--r')
ax5.set_xlabel('genomic_position', fontdict=ax_fontdict)
ax5.set_ylabel('nucleotide diversity', fontdict=ax_fontdict)
ax6 = fig.add_subplot(246)
# plt.title(('Nucleotide diversity 100 timesteps into sweep\n(calculated in '
#          '11-locus windows'))
ax6.plot(range(mod.comm[0].gen_arch.L), pi[1])
ax6.plot([nonneut_loc, nonneut_loc], [0, 1], '--r')
ax6.set_xlabel('genomic_position', fontdict=ax_fontdict)
# ax6.set_ylabel('nucleotide diversity', fontdict=ax_fontdict)
ax6.set_yticks([])
ax7 = fig.add_subplot(247)
# plt.title(('Nucleotide diversity at end of sweep\n(calculated in 11-locus '
#          'windows)'))
ax7.plot(range(mod.comm[0].gen_arch.L), pi[2])
ax7.plot([nonneut_loc, nonneut_loc], [0, 1], '--r')
ax7.set_xlabel('genomic_position', fontdict=ax_fontdict)
ax7.set_yticks([])
# ax7.set_ylabel('nucleotide diversity', fontdict=ax_fontdict)
ax8 = fig.add_subplot(248)
# plt.title(('Nucleotide diversity 2500 timesteps after end of sweep\n'
#           '(calculated in 11-locus windows)'))
ax8.plot(range(mod.comm[0].gen_arch.L), pi[3])
ax8.plot([nonneut_loc, nonneut_loc], [0, 1], '--r')
ax8.set_xlabel('genomic_position', fontdict=ax_fontdict)
ax8.set_yticks([])
# ax8.set_ylabel('nucleotide diversity', fontdict=ax_fontdict)
pi_min_lim = 0.95 * min([val for sublist in pi for val in sublist])
pi_max_lim = 1.05 * max([val for sublist in pi for val in sublist])
ax5.set_ylim((pi_min_lim, pi_max_lim))
ax6.set_ylim((pi_min_lim, pi_max_lim))
ax7.set_ylim((pi_min_lim, pi_max_lim))
ax8.set_ylim((pi_min_lim, pi_max_lim))

# add vertical space between the first and second rows of plots
fig.subplots_adjust(hspace=.5)

plt.show()
plt.savefig(os.path.join(img_dir, 'SWEEP_pop_and_nuc_div.pdf'))
print('Sweep complete')

plt.rcParams['figure.figsize'] = [5.5, 4]
fig2 = plt.figure()
ax = fig2.add_subplot(111)
# plt.title('Mean fitness across model time', fontdict=ttl_fontdict)
ax.plot(np.linspace(0, mod.t, len(mean_fit)), mean_fit, '-', color='#D55E00')
ax.set_xlabel('time', fontdict=ax_fontdict)
ax.set_ylabel('mean fitness', fontdict=ax_fontdict)
ax.set_ylim(((1 - mod.comm[0].gen_arch.traits[0].phi)-0.05, 1.01))
ax.set_xlim((0, mod.t))
ax.set_yticks(np.linspace(0.8, 1.0, 5))

plt.show()
plt.savefig(os.path.join(img_dir, 'SWEEP_mean_fit.pdf'))
