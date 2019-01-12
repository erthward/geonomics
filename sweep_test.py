#!/usr/bin/python
# sweep_test.py

import geonomics as gnx

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter as C
from itertools import combinations


def calc_pi(mod, window_width=5):
    print("\nCalculating nucleotide diversity (WILL TAKE A FEW MINUTES)\n")
    # make speciome
    speciome = np.stack([i.genome for i in mod.comm[0].values()])
    # create data structure to store pi value for each genomic window
    pi_vals = []
    # for each locus
    for loc in range(speciome.shape[1]):
        # get 5-locus window around the locus
        start = max(loc-window_width, 0)
        # add 6 so that the resulting interval includes 5 loci on either
        # side of focal locus
        stop = min(loc+window_width+1, speciome.shape[1]-1)
        print('\tcalculating for interval [%i, %i)' % (start, stop))
        # get all individs' genotypes for that interval of loci
        intervalome = speciome[:, range(start, stop), :]
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
                [i.genome[nonneut_loc, :] for i in mod.comm[0].values()]))
    return nonneut_loc_freq


# create data structures to save data to be analyzed/plotted
mean_fit = []
pi = []

# build the model
mod = gnx.make_model('./test/validation/sweep/sweep_params.py')
# set the non-neutral locus to the central locus (i.e. 50, in a 101-length
# genomic arch.), then set the starting 1-allele freq of the
# non-neutral allele to 0
nonneut_loc = 50
mod.comm[0].gen_arch.traits[0].loci[0] = nonneut_loc
mod.comm[0].gen_arch.p[nonneut_loc] = 0
# burn the model in, then run for 50 timesteps (to make sure pop size is
# at equilibrium)
mod.walk(20000, 'burn', True)
mod.walk(50, 'main', True)

# create the figure and first Axes instance
fig = plt.figure()
ax1 = fig.add_subplot(231)
ax1.set_title('Population at beginning of sweep\n(colored by phenotype)')

# keep reintroducing a mutant and running for 100 timesteps, until
# mutant establishes
sweeping_allele_established = False
lapsed_t = 0
while not sweeping_allele_established:

    mean_fit = [-999]
    pi = []

    # store the starting mean fitness
    mean_fit[0] = calc_mean_fit(mod)

    # introduce the adaptive allele in a random individual
    print('\n=================\nINTRODUCING MUTATION\n')
    individ = np.random.choice([*mod.comm[0]])
    chromatid = np.random.binomial(1, 0.5, 1)
    mod.comm[0][individ].genome[nonneut_loc, chromatid] = 1

    # plot the beginning of the sweep, with the mutant individual
    plt.gca()
    mod.plot_phenotype(0, 0, 0, size=10)
    ax1.plot(mod.comm[0][individ].x - 0.5, mod.comm[0][individ].y - 0.5, 'ow')

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
    else:
        print('\n\nMutation lost...\n')

    # increment lapsed_t
    lapsed_t += 50

# calculate and store nucleotide diversity 50 timesteps into the sweep
pi.append(calc_pi(mod))

# plot status after 100 more timesteps
for _ in range(10):
    mod.walk(10)
    # calclate and store mean fitness
    mean_fit.append(calc_mean_fit(mod))
    print('\tallele_freq = %0.3f' % calc_nonneut_loc_freq(mod))
ax2 = fig.add_subplot(232)
plt.title('Population 200 timesteps into sweep')
mod.plot_phenotype(0, 0, 0, size=10)

# walk until mutant fixed
nonneut_loc_freq = np.mean(np.vstack([i.genome[nonneut_loc,
                                               :] for i in mod.comm[
                                                                0].values()]))
fixed = nonneut_loc_freq == 1
while not fixed:
    mod.walk(10, verbose=True)
    print('\tallele_freq = %0.3f' % calc_nonneut_loc_freq(mod))
    # calclate and store mean fitness
    mean_fit.append(calc_mean_fit(mod))
    # check if fixed
    nonneut_loc_freq = np.mean(np.vstack(
                               [i.genome[nonneut_loc,
                                         :] for i in mod.comm[0].values()]))
    fixed = nonneut_loc_freq == 1

# plot status after sweep is complete
ax3 = fig.add_subplot(233)
plt.title('Population after sweep has completed (t = %i)' % (
                                                mod.t - lapsed_t + 50 + 50))
mod.plot_phenotype(0, 0, 0, size=10)

# calculate and store ending nucleotide diversity
pi.append(calc_pi(mod))

# plot stats
ax4 = fig.add_subplot(234)
plt.title('Mean fitness across model time')
ax4.plot(range(len(mean_fit)), mean_fit)
ax4.set_xlabel('model time')
ax4.set_ylabel('mean fitness')
ax4.set_ylim((1 - mod.comm[0].gen_arch.traits[0].phi, 1.0))
ax5 = fig.add_subplot(235)
plt.title(('Nucleotide diversity 100 timesteps into sweep\n(calculated in '
          '11-locus '))
ax5.plot(range(mod.comm[0].gen_arch.L), pi[0])
ax5.plot([50, 50], [0, 1], '--r')
ax5.set_xlabel('genomic_position')
ax5.set_ylabel('$\\Pi$')
ax6 = fig.add_subplot(236)
plt.title(('Nucleotide diversity at end of sweep\n(calculated in 11-locus '
          'windows)'))
ax6.plot(range(mod.comm[0].gen_arch.L), pi[1])
ax6.plot([50, 50], [0, 1], '--r')
ax6.set_xlabel('genomic_position')
ax6.set_ylabel('$\\Pi$')
pi_min_lim = 0.95 * min([val for sublist in pi for val in sublist])
pi_max_lim = 1.05 * max([val for sublist in pi for val in sublist])
ax5.set_ylim((pi_min_lim, pi_max_lim))
ax6.set_ylim((pi_min_lim, pi_max_lim))
plt.show()
print('Sweep complete')
