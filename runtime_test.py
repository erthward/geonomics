#!usr/bin/python

# runtime_test.py

import geonomics as gnx
import numpy as np
from copy import deepcopy
import time
import matplotlib.pyplot as plt
import os

# set some image params
img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
           'img/final/')
ax_fontdict = {'fontsize': 12,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': 16,
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
             'n_births_distr_lambda': 2}

# define the ranges of values for the parameters I want to vary
test_vals = {'L': [10, 100, 1000],
             'K_fact': [1, 5, 10, 20],
             'dim': [(10, 10), (20, 20), (50, 50), (100, 100)],
             'n_births_distr_lambda': [0, 1, 2],
             }

# define a data structure to save my runtimes in
runtimes = {}
# and a data structure to save mean population sizes
mean_pop_sizes = {}
std_pop_sizes = {}

# start my figure
plt.rcParams['figure.figsize'] = [6, 2]
fig = plt.figure()
fig.suptitle(('Geonomics runtime analysis\n(averaged over %i '
              'timesteps per parameterization)') % T)

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
        mod.walk(1000000, 'burn', verbose=False)
        print('\tSTARTING MODEL...')
        start = time.time()
        # run the model
        mod.walk(T, 'main', verbose=False)
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

    # plot my runtimes
    if n == 0:
        ax = fig.add_subplot(141)
        ax0 = ax
    else:
        ax = fig.add_subplot(141 + n, sharey=ax0)
    ax.set_ylim((0, 60))
    # if param == 'K_fact':
    #    plt.plot(mean_pop_sizes_this_param, runtimes_this_param, '-b')
    #    ax.set_xlabel('mean population size (as a result of K_factor)',
    #                  fontdict=ax_fontdict)
    # elif param == 'dim':
    if param == 'dim':
        plt.plot([i[0] for i in test_vals[param]], runtimes_this_param, '-or')
    else:
        plt.plot(test_vals[param], runtimes_this_param, '-or')
    ax.set_xlabel(param, fontdict=ax_fontdict)
    ax.set_ylabel('mean runtime (sec/timestep)', fontdict=ax_fontdict)

    # save all the runtimes for this param in the runtimes dict
    runtimes[param] = runtimes_this_param
    mean_pop_sizes[param] = mean_pop_sizes_this_param
    std_pop_sizes[param] = std_pop_sizes_this_param

plt.show()
plt.savefig(os.path.join(img_dir, 'runtime.pdf'))
