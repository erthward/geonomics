#!/usr/bin/python

# flake8: noqa


# bottleneck_test.py

import geonomics as gnx

import os
import re
import shutil
import numpy as np
import pandas as pd
import vcf
import matplotlib.pyplot as plt

# set some plotting params
#img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
           'img/final/')
titlesize=20
axlabelsize=18
ticklabelsize=15

# set the data directory, and delete it if it already exists (so that we don't
# create mutliple, conflicting sets of data files)
data_dir = './GNX_mod-bottleneck_params'
if os.path.isdir(data_dir):
    shutil.rmtree(data_dir)

# make and run model
mod = gnx.make_model('./tests/validation/bottleneck/bottleneck_params.py')
mod.run(verbose=True)

# for each iteration
its_dirs = os.listdir(data_dir)
for it_dir in its_dirs:
    # read in the data
    files_dir = os.path.join(data_dir, it_dir, 'spp-spp_0')
    files = os.listdir(files_dir)
    vcf_files = [f for f in files if os.path.splitext(f)[1] == '.vcf']
    csv_files = [f for f in files if re.search('spp_0\\.csv$', f)]
    timesteps = [re.search('(?<=t\\-)\\d*(?=\\_spp)', f).group(
                    ) for f in vcf_files]
    timesteps = sorted([int(step) for step in timesteps])
    # create data structure to store allele frequencies
    allele_freqs = {loc: [] for loc in range(mod.comm[0].gen_arch.L)}
    for step in timesteps:
        vcf_file = [f for f in vcf_files if re.search(
                            '(?<=t\\-)%i(?=\\_spp)' % step, f)]
        assert len(vcf_file) == 1, 'More than one VCF file!'
        vcf_file = vcf_file[0]
        vcf_reader = vcf.Reader(open(os.path.join(files_dir, vcf_file), 'r'))
        csv_file = [f for f in csv_files if re.search(
                            '(?<=t\\-)%i(?=\\_spp)' % step, f)]
        assert len(csv_file) == 1, 'More than one CSV file!'
        csv_file = csv_file[0]
        csv = pd.read_csv(os.path.join(files_dir, csv_file))
        # grab all individuals' genotypes into an n_individs x n_loci array
        genotypes = np.ones((len(csv), mod.comm[0].gen_arch.L))*99
        try:
            for loc in range(genotypes.shape[1]):
                rec = next(vcf_reader)
                for n, ind in enumerate(csv['idx']):
                    genotype = sum([int(base) for base in rec.genotype(
                        str(ind)).data.GT.split('|')])/2
                    genotypes[n, loc] = genotype
        except StopIteration:
            pass
        assert 99 not in genotypes
        freqs = np.sum(genotypes, axis=0) / (genotypes.shape[0])
        for loc in range(len(freqs)):
            allele_freqs[loc].append(freqs[loc])


# plot results
spp = mod.comm[0]
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
fig = plt.figure()
#fig.suptitle('Drift during bottleneck event', size=20)
ax1 = fig.add_subplot(311)
ax1.set_title('allele-frequency trajectories', size=titlesize)
ax1.set_xlabel('time', size=axlabelsize)
ax1.set_ylabel("frequency of '1' allele", size=axlabelsize)
ax1.tick_params(labelsize=ticklabelsize)
plt.xlim((0, mod.T))
plt.ylim((0, 1))
for loc, traj in [*allele_freqs.items()][:10]:
    plt.plot(timesteps, traj, '-k', linewidth=1, alpha=0.5)
# add lines for start and end of the bottleneck event
# (dem-change events take place as changes in the carrying-capacity rasters
# at the ends of their scheduled timesteps, after pop dynamics have happened,
# (in this case timestpes 49 and 54) so makes most sense to plot
# the lines at the following timesteps (50 and 55),
# which are the first timeteps in which the events have an effect
# on the number of individuals)
ts = mod.params.comm.species['spp_0'].change.dem[0].timesteps
plt.plot([ts[0]] * 2, [0, max(spp.Nt)], '--r')
plt.plot([ts[1]] * 2, [0, max(spp.Nt)], '--r')
ax2 = fig.add_subplot(312)
ax2.set_title('population size', size=titlesize)
ax2.set_xlabel('time', size=axlabelsize)
ax2.set_ylabel('population size', size=axlabelsize)
ax2.tick_params(labelsize=ticklabelsize)
plt.xlim((0, mod.T))
plt.ylim((0, 1.1*max(spp.Nt)))
plt.plot(range(mod.T), spp.Nt[-mod.T:], '-k')
# add lines for start and end of the bottleneck event
plt.plot([ts[0]] * 2, [0, 1.15 * max(spp.Nt)], '--r')
plt.plot([ts[1]] * 2, [0, 1.15 * max(spp.Nt)], '--r')
plt.show()


# estimate effective population size using harmonic and arithmetic means
Ne_harm = 1 / ((1 / len(spp.Nt)) * sum([1 / N for N in spp.Nt]))
Ne_arit = np.mean(spp.Nt)
# estimate mean times to fixation for pops with each N_e, for loci starting at
# p = q = 0.5 (according to W-F model)
t_fix_harm = 2.776 * Ne_harm
t_fix_arit = 2.776 * Ne_arit
# TODO: compare both expected times to fixation to the mean time
# to fixation in the simulation


# calculate allele-freq change rates in moving windows
window_width = 3
window_buff = int(window_width / 2)
change_rates = {}
for i in range(window_buff,
               max([len(freq_list) for freq_list in allele_freqs.values(
                                                            )]) - window_buff):
    loc_change_rates = []
    for loc in range(len(allele_freqs)):
        loc_change = (allele_freqs[loc][i + window_buff] -
                      allele_freqs[loc][i - window_buff])
        loc_change_rate = np.abs(loc_change) / window_width
        loc_change_rates.append(loc_change_rate)
    change_rate = np.mean(loc_change_rates)
    change_rates[timesteps[i]] = change_rate

ax3 = fig.add_subplot(313)
ax3.set_title('allele-frequency change rates', size=titlesize)
ax3.set_xlabel('time', size=axlabelsize)
ax3.set_ylabel('rate of allele-frequency change', size=axlabelsize)
ax3.tick_params(labelsize=ticklabelsize)
plt.plot([*change_rates.keys()], [*change_rates.values()], '-k')
plt.plot([ts[0]] * 2, [0, 1.1 * max([*change_rates.values()])], '--r')
plt.plot([ts[1]] * 2, [0, 1.1 * max([*change_rates.values()])], '--r')
plt.xlim((0, 300))
plt.ylim((0, 1.05 * max([*change_rates.values()])))
plt.show()

#plt.savefig(os.path.join(img_dir, 'BOTTLENECK_bottleneck_plot.pdf'))
