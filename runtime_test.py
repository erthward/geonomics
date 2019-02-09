#!usr/bin/python

# runtime_test.py

import geonomics as gnx
import numpy as np
from copy import deepcopy
import time
import matplotlib.pyplot as plt

# read in the parameters file
params = gnx.read_parameters_file(('./geonomics/tests/runtime/'
                                   'runtime_params.py'))

# define the base values of the parameters (i.e. the values at which they'll
# be held in all simulations where they are not the parameter whose influence
# on runtime is being tested)
base_vals = {'K_fact': 1,
             'dim': (50, 50),
             'L': 1000,
             'n_births_distr_lambda': 2}

# define the ranges of values for the parameters I want to vary
test_vals = {'K_fact': [1, 5, 10, 100],
             'dim': [(10, 10), (20, 20), (50, 50), (100, 100)],
             'L': [10, 100, 1000, 10000],
             'n_births_distr_lambda': [0, 1, 2, 5],
             }

# define a data structure to save my runtimes in
runtimes = {}

# start my figure
fig = plt.figure()
fig.suptitle('Geonomics runtime analysis')

# for each parameter I want to test various values of
for n, param in enumerate([*test_vals]):
    # get the names of the other parameters, to be held constant
    other_params_list = [p for p in [*test_vals] if p is not param]

    # create a list to save the runtimes for this param
    runtimes_this_param = []

    # get a deepcopy of the params object for the runtime analyses for this
    # param
    new_params = deepcopy(params)

    # define a dictionary of lambda functions to reset each parameter's values
    # in the new_params object
    reset_fns = {'dim': lambda dim: setattr(
                    new_params['landscape']['main'],
                 'dim', dim),
                 'K_fact': lambda K_fact: setattr(
                    new_params['comm']['species']['spp_0']['init'],
                    'K_factor', K_fact),
                 'L': lambda L: setattr(
                    new_params['comm']['species']['spp_0']['init'],
                    'K_factor', L),
                 'n_births_distr_lambda': lambda lambda_val: setattr(
                    new_params['comm']['species']['spp_0']['mating'],
                    'n_births_distr_lambda', lambda_val)

                 }

    # set the 'other params' (the ones whose effect on runtime are not being
    # tested in this loop) to their base values
    for other_param in other_params_list:
        reset_fns[other_param](base_vals[other_param])

    # for each value of the param being tested
    for param_val in test_vals[param]:

        # reset the parameter's value in the new_params object
        reset_fns[param](param_val)

        # reset the layer_0 rast, if param is dim
        if param == 'dim':
            lyr_params = new_params['landscape']['layers']
            setattr(lyr_params['layer_0']['init']['defined'], 'rast',
                    np.ones(param_val))
        # reset the chrom-length list if param is L
        if param == 'L':
            gen_params = new_params['comm']['species']['spp_0']['gen_arch']
            setattr(gen_params, 'l_c', [param_val])

        # get the start time
        start = time.time()

        # make the model
        mod = gnx.make_model(new_params)

        # run the model
        mod.run(verbose=False)

        # get the stop time, calculate the elapsed time, and append it to list
        stop = time.time()
        elapsed = stop-start
        runtimes_this_param.append(elapsed)

        # plot my runtimes
        ax = fig.addsubplot(141 + n)
        plt.plot(test_vals[param], runtimes_this_param)
        ax.set_xlabel(param)
        ax.set_ylabel('runtime (seconds)')

    # save all the runtimes for this param in the runtimes dict
    runtimes[param] = runtimes_this_param
