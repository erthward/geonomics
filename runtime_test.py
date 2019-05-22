#!usr/bin/python

# runtime_test.py

import geonomics as gnx
import numpy as np
from copy import deepcopy
import time
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

# set some image params
img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
           'img/final/')
ax_fontdict = {'fontsize': 10,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': 12,
                'name': 'Bitstream Vera Sans'}

# read in the parameters file
params = gnx.read_parameters_file(('./geonomics/tests/runtime/'
                                   'runtime_params.py'))

# define number of timesteps to run for each model
T = 250

# define the base values of the parameters (i.e. the values at which they'll
# be held in all simulations where they are not the parameter whose influence
# on runtime is being tested)
base_vals = {'L': 1000,
             'K_fact': 5,
             'dim': (20, 20),
             'n_births_distr_lambda': 2
             }
# define the ranges of values for the parameters I want to vary
test_vals = {'L': [10, 100, 1000],
             'K_fact': [5, 10, 20],
             'dim': [(10, 10), (20, 20), (50, 50), (100, 100)],
             'n_births_distr_lambda': [0, 1, 2],
             }
# define x-axis labels for plots
plot_x_labs = {'L': 'genome\nlength (L)',
               'K_fact': 'carrying-capacity\nconstant (K)',
               'dim': 'landscape\nsize (dim)',
               'n_births_distr_lambda': ('mean n. births per\n '
                                         'mating pair (lambda)')
               }

# define a data structure to save my runtimes in
runtimes = {}
# and a data structure to save mean population sizes
mean_pop_sizes = {}
std_pop_sizes = {}

# start my figure
fig = plt.figure(figsize=(6.75, 3))
gs = gridspec.GridSpec(1, 4)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])
ax4 = plt.subplot(gs[3])
axs = [ax1, ax2, ax3, ax4]
#plt.rc('xtick',labelsize=8)
#plt.rc('ytick',labelsize=8)
#plt.rcParams['xtick.labelsize']=2
#plt.rcParams['ytick.labelsize']=2
# fig.suptitle(('Geonomics runtime analysis\n(averaged over %i '
#              'timesteps per parameterization)') % T)

# for each parameter I want to test various values of
for n, param in enumerate([*test_vals]):
    print('Now testing runtimes across various values of %s...\n' % param)
    # get the names of the other parameters, to be held constant
    other_params_list = [p for p in [*test_vals] if p is not param]

    # create a list to save the runtimes and mean pop sizes for this param
    runtimes_this_param = []
    mean_pop_sizes_this_param = []
    std_pop_sizes_this_param = []

    # get a deepcopy of the params object for the runtime analyses for this
    # param
    new_params = deepcopy(params)

    # define a dictionary of lambda functions to reset each parameter's values
    # in the new_params object
    reset_fns = {'dim': lambda dim: new_params['landscape'][
                                'main'].update({'dim': dim}),
                 'K_fact': lambda K_fact: new_params['comm']['species'][
                                'spp_0']['init'].update({'K_factor': K_fact}),
                 'L': lambda L: new_params['comm']['species']['spp_0'][
                                'gen_arch'].update({'L': L}),
                 'n_births_distr_lambda': lambda lambda_val: new_params[
                                'comm']['species']['spp_0'][
                                'mating'].update(
                                        {'n_births_distr_lambda': lambda_val})
                 }

    # set the 'other params' (the ones whose effect on runtime are not being
    # tested in this loop) to their base values
    for other_param in other_params_list:
        reset_fns[other_param](base_vals[other_param])

    # for each value of the param being tested
    for param_val in test_vals[param]:
        print('\tTesting value %s...' % str(param_val))

        # reset the parameter's value in the new_params object
        reset_fns[param](param_val)

        # reset the layer_0 rast, if param is dim
        if param == 'dim':
            new_params['landscape']['layers']['layer_0']['init']['defined'][
                                    'rast'] = np.ones(param_val)

        # reset the chrom-length list if param is L
        if param == 'L':
            new_params['comm']['species']['spp_0']['gen_arch'][
                                    'l_c'] = [param_val]

        print('\t' + '-'*70)
        print('\tK_factor:')
        print('\t\t' + str(new_params.comm.species.spp_0.init.K_factor))
        print('\tdim')
        print('\t\t' + str(new_params.landscape.main.dim))
        print('\tL')
        print('\t\t' + str(new_params.comm.species.spp_0.gen_arch.L))
        print('\tn_births_distr_lambda')
        print('\t\t' + str(
            new_params.comm.species.spp_0.mating.n_births_distr_lambda))
        print('\t' + '-'*70)

        # get the start time
        # make the model
        mod = gnx.make_model(new_params)
        print('\tMODEL MADE')
        # burn it in
        mod.walk(1000000, 'burn', verbose=True)
        print('\tSTARTING MODEL...')
        start = time.time()
        # run the model
        mod.walk(T, 'main', verbose=True)
        # get the stop time, calculate the elapsed time, and append it to list
        stop = time.time()
        elapsed = stop-start
        mean_runtime = elapsed/T
        print('\tMODEL FINISHED')
        print('\tmean runtime: %0.6f' % mean_runtime)
        runtimes_this_param.append(mean_runtime)
        # calculate and append the mean and std population size
        mean_pop_size = np.mean(mod.comm[0].Nt)
        print('\tmean pop size: %0.2f' % mean_pop_size)
        mean_pop_sizes_this_param.append(mean_pop_size)
        std_pop_size = np.std(mod.comm[0].Nt)
        print('\tstd pop size: %0.2f\n\n\n\n' % std_pop_size)
        std_pop_sizes_this_param.append(std_pop_size)

    if param == 'dim':
        axs[n].plot([i[0] for i in test_vals[param]],
                    runtimes_this_param, '-or')
    else:
        axs[n].plot(test_vals[param], runtimes_this_param, '-or')
    axs[n].set_xlabel(plot_x_labs[param], fontdict=ax_fontdict)

    # save all the runtimes for this param in the runtimes dict
    runtimes[param] = runtimes_this_param
    mean_pop_sizes[param] = mean_pop_sizes_this_param
    std_pop_sizes[param] = std_pop_sizes_this_param



assert True == False, ("BECAUSE THE PLOT DOESN'T FORMAT CORRECTLY WHEN THE "
                       "SCRIPT IS RUN ALL THE WAY ON ITS OWN, IT BREAKS HERE "
                       "SO THAT YOU CAN RUN THE REMAINING LINES MANUALLY TO "
                       "FORMAT AND SAVE THE PLOT.")
# some final plot tweaking
ax1.set_ylabel('mean runtime\n(sec/timestep)',
               fontdict=ax_fontdict)
ymax = 1.1 * max([max(rts) for rts in runtimes.values()])
ax1.set_ylim((0, ymax))
ax2.set_ylim((0, ymax))
ax3.set_ylim((0, ymax))
ax4.set_ylim((0, ymax))
ax2.set_yticks([])
ax3.set_yticks([])
ax4.set_yticks([])
ax1.tick_params(size=3)
ax2.tick_params(size=3)
ax3.tick_params(size=3)
ax4.tick_params(size=3)
ax1.set_xticklabels(ax1.get_xticklabels(), fontdict={'size': 8})
ax2.set_xticklabels(ax2.get_xticklabels(), fontdict={'size': 8})
ax3.set_xticklabels(ax3.get_xticklabels(), fontdict={'size': 8})
ax4.set_xticklabels(ax4.get_xticklabels(), fontdict={'size': 8})
ax1.set_yticklabels(ax1.get_yticklabels(), fontdict={'size': 8})
ax2.set_yticklabels(ax2.get_yticklabels(), fontdict={'size': 8})
ax3.set_yticklabels(ax3.get_yticklabels(), fontdict={'size': 8})
ax4.set_yticklabels(ax4.get_yticklabels(), fontdict={'size': 8})
plt.subplots_adjust(bottom=0.3,
                    top=0.88,
                    left=0.12,
                    right=0.90,
                    wspace=0.20,
                    hspace=0.20)
plt.show()
plt.savefig(os.path.join(img_dir, 'runtime.pdf'),
            format='pdf',
            dpi=1000)
