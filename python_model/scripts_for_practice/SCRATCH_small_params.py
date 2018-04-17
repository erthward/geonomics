#!/usr/bin/python

import numpy as np

params = {
    'set_seed': True,
    'seed_num': 2,
    'T': 2500,  # total model runtime
    'burn_T_min': 50,  # total burn-in runtime
    'L': 1000,
    'l_c': [500, 250, 250],
    'recomb_array': None,
    'x': 2,
    'mu': 10e-9,
    'alpha_r': 0.5,
    'beta_r': 400,
    'alpha_mut_s': 25,
    'beta_mut_s': 0.5,
    'use_dom': False,
    'pleiotropy': True,
    'traits': {'num': 2,
               'scape_num': [0, 0],
               'phi': [0.05, 0.05],
               'n': [1, 100],
               'alpha_dist_sigma': [0.5, 0.00025],
               'fitness_fn_gamma': [1, 1],
               'univ_advant': [False, False]
               },

    'N': 10,  # total pop size

    'dims': (3, 3),  # dimensions of landscape

    'density_grid_window_width': None,  # with window-width used for the Density_Grid_Stack used to calculate
    # population density (if set to None, defaults to the closest factor
    # of the larger landscape dimension to 1/10th of that dimension)

    'num_scapes': 2,  # number of landscapes desired

    'rand_land': True,

    'n_rand_pts': 2600,
    # number of random coordinates to be used in generating random landscapes (only needed if rand_land = True)

    'islands': True,

    'island_val': 0.1,
    'island_mask': None,
    'landscape_pt_coords': np.array(
        [[0, 0], [0, 100], [100, 0], [50, 40], [100, 100], [30, 0], [0, 30], [70, 100], [100, 70]]),
    # coords of points to use to interpolate defined landscape layers (can be provided as either a single nx2 Numpy
    # array, where n matches the number of points in landscape_pt_vals arrays, to be used as the points for each
    # landscape layer, or a list or tuple of nx2 Numpy arrays, one for each landscape layer; only needed if rand_land
    #  = False)

    'landscape_pt_vals': [np.array([-0.5, 0.5, 0.5, 0.5, 1.5, -0.5, -0.5, 1.5, 1.5]),
                          np.array([1.5, -0.5, -0.5, -0.5, 1.5, 0.85, 0.85, 0.85, 0.85])],
    # point values to use to interpolate defined landscape layers (a list or tuple of 1xn Numpy arrays,
    # where n matches the number of points in landscape_pt_coords arrays; only needed in rand_land = False)

    'interp_method': ['linear', 'linear'],

    'K_cap': 0.4,  # per-cell highest carrying capacity value to be reached during burn-in

    'move': True,  # is this a mobile species?

    'movement_surf': True,
    # use a landscape layer as a resistance surface, or habitat quality layer, to direct movement?
    # 'movement_surf' : False,

    'n_movement_surf_scape': 1,

    'movement_surf_vonmises_kappa': 2,

    'movement_surf_gauss_KDE_bandwidth': 0.2,

    'mu_direction': 0,

    'kappa_direction': 0,

    'mu_distance': 0.5,

    'sigma_distance': 0.5,

    'sex': False,

    'repro_age': 0,

    'dist_weighted_birth': False,

    'R': 0.5,  # pop intrinsic growth rate

    'b': 0.2,
    # population intrinsic birth rate (implemented as the probability that an identified potential mating pair
    # successfully mates); NOTE: this may later need to be re-implemented to allow for spatial variation in intrinsic
    #  rate (e.g. expression as a raster) and/or for density-dependent births as well as deaths

    'd_min': 0.01,  # minimum neutral (i.e. non-selection driven) probability of death

    'd_max': 0.90,  # maximum neutral probability of death

    'lambda_offspring': 4,
    # expected value of offspring for a successful mating pair (used as the lambda value in a Poisson distribution)

    'mating_radius': 0.5,  # radius of mate-searching area

    'mu_dispersal': 0.5,  # mean dispersal distance (lognormal distribution)

    'sigma_dispersal': 0.5,  # sd of dispersal distance

    'size': 1,
    # float/int, or list/tuple of length T containing floats/ints, expressing the target population size over model
    # time as a ratio of the starting size (N)

    'sigma_deaths': 0.2,
    # std for the normal distribution used to choose the r.v. of deaths per timestep (mean of this distribution is the
    # overshoot, as calculated from pop.size and pop.census())

    'density_dependent_fitness': True,
    # should fitness be density dependent? (note: helps to avoid subpopulation 'clumping')

    'custom_fns': {'recomb_rate_fn': None
                   # if provided, must be a function that returns a single recombination rate value (r) when called

                   },

    'data': {
        # dictionary defining the data to be collected, the sampling strategy to use, the timesteps for collection,
        # and other parameters

        'sampling_scheme': 'random',  # can be 'all', 'random', 'point', or 'transect'
        'sampling_args': {'n': 50
                          },

        # args to be unpacked into sampling function (see docstring of sample_data function in data module for details)
        'freq': 15,
        # can be an integer (in which case data will be collected every that many timesteps, plus at the end)
        # or a list of specific timesteps

        'include_land': False,
        # if True, will save the Landscape_Stack object each time other data is saved (probably only useful if land
        # is changing in some way not manually coded by the user)

        'gen_data_format': 'VCF',  # can be 'VCF', 'FASTA', or 'ms'

        'geo_data_format': ['CSV', 'Geotiff'],
        # 1st argument for points, 2nd for raster;
        # currently 1.) CSV, Shapefile and 2.) Geotiff available

        'run_name': 'test',  # a name for this parameterization and run (used to name the data_directory and files)
        'write_intermittent': True,
        'drop_after_write': True
    },

    'stats': {
        # dictionary defining which stats to be calculated, and parameters on their calculation (including frequency,
        #  in timesteps, of collection) valid stats include: 'Nt'  : population census size 'het' : heterozygosity
        # 'MAF' : minor allele frequency 'LD'  : linkage disequilibrium
        'Nt': {'calc': True, 'freq': 1},

        'het': {'calc': True, 'freq': 1},

        'MAF': {'calc': True, 'freq': 1},

        'LD': {'calc': True, 'freq': 1}

    }

}
