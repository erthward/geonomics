#!/usr/bin/python
# sweep_test.py

import geonomics as gnx

import numpy as np
import matplotlib.pyplot as plt


# define function to plot a timestep
def plot_timestep(mod, t, mean_fit):
    fig = plt.figure()
    plt.suptitle('Selective sweep, timestep=%i' % t)
    ax1 = fig.add_subplot(131)
    ax1.set_title('Mean fitness over time')
    plt.xlabel('time')
    plt.ylabel('population mean fitness')
    ax1.plot(range(len(mean_fit)), mean_fit, '-k')
    ax3 = fig.add_subplot(133)
    ax3.set_title('Major allele frequencies for all loci')
    plt.xlabel('genomic position')
    plt.ylabel('major allele frequency')
    mafs = []
    for loc in range(mod.comm[0].gen_arch.L):
        speciome = np.vstack([i.genome[loc, :] for i in mod.comm[0].values()])
        freq_1_allele = np.mean(speciome)
        maf = max(freq_1_allele, 1 - freq_1_allele)
        mafs.append(maf)
    ax3.plot(range(len(mafs)), mafs)
    plt.ylim((0, 1))
    ax2 = fig.add_subplot(132)
    ax2.set_title(('Population colored by phenotype (outer circle)\n'
                   'and fitness (inner circle)'))
    mod.plot_fitness(0, 0, 0, fitness_cbar=False)
    plt.savefig('sweep_plot_timestep_%i' % t)
    plt.close()
    return mafs


def calc_mean_fit(mod):
    mean_fit = np.mean(mod.comm[0]._get_fit())
    return mean_fit


def run_mod_n_steps(t, n, mod, mean_fit, mafs):
    for t in range(n):
        mod.walk(1, 'main', True)
        mean_fit.append(calc_mean_fit(mod))
    t += n
    mafs_this_timestep = plot_timestep(mod=mod, t=t, mean_fit=mean_fit)
    mafs.append(mafs_this_timestep)
    return(t, mean_fit, mafs)


# run the whole thing for 300 timesteps, then abort and restart
# if the introduced allele is lost (otherwise will continue to run
# the model until the mutation has fixed)
sweeping_allele_established = False

while not sweeping_allele_established:

    # create data structures to save data to be analyzed/plotted
    mean_fit = []
    mafs = []
    t = 0

    # build the model
    mod = gnx.make_model('./test/validation/sweep/sweep_params.py')

    # set the non-neutral locus to the central locus (i.e. 50, in a 101-length
    # genomic arch.), then set the starting 1-allele freq of the
    # non-neutral allele to 0
    nonneut_loc = 50
    mod.comm[0].gen_arch.traits[0].loci[0] = nonneut_loc
    mod.comm[0].gen_arch.p[nonneut_loc] = 0

    # burn the model in
    mod.walk(20000, 'burn', True)

    # run the model for 1 main timestep (to have some data to plot in the first
    # plot)
    mod.walk(1, 'main', True)
    t += 1

    # store the mean fitness and MAFs after this timestep
    mean_fit.append(calc_mean_fit(mod))

    # plot population and data after the first timestep
    mafs_this_timestep = plot_timestep(mod=mod, t=t, mean_fit=mean_fit)
    mafs.append(mafs_this_timestep)

    # walk it for 100 timesteps, plotting mean fitness the whole time
    t, mean_fit, mafs = run_mod_n_steps(t, 100, mod, mean_fit, mafs)

    # introduce the adaptive allele in a random individual
    individ = np.random.choice([*mod.comm[0]])
    chromatid = np.random.binomial(1, 0.5, 1)
    mod.comm[0][individ].genome[nonneut_loc, chromatid] = 1

    # walk it for 200 more timesteps, plotting mean fitness the whole time
    t, mean_fit, mafs = run_mod_n_steps(t, 100, mod, mean_fit, mafs)

    # check if the allele was established
    nonneut_loc_sum = sum(np.vstack(
                [i.genome[nonneut_loc, :] for i in mod.comm[0].values()]))
    sweeping_allele_established = nonneut_loc_sum > 0
    print(('\n\nSWEEPING ALLELE DID%s ESTABLISH.'
          '\n\n') % (' NOT' * (not sweeping_allele_established)))

nonneut_loc_freq = np.mean(np.vstack(
                           [i.genome[nonneut_loc,
                                     :] for i in mod.comm[0].values()]))
fixed = nonneut_loc_freq == 1
print("AFTER ESTABLISHING ALLELE, t = %i" % t)
while not fixed:
    mod.walk(1, verbose=True)
    nonneut_loc_freq = np.mean(np.vstack(
                               [i.genome[nonneut_loc,
                                         :] for i in mod.comm[0].values()]))
    fixed = nonneut_loc_freq == 1
    mean_fit.append(calc_mean_fit(mod))
    if t % 100 == 0:
        mafs_this_timestep = plot_timestep(mod=mod, t=t, mean_fit=mean_fit)
        mafs.append(mafs_this_timestep)
    t += 1

print('Sweep occurred')
