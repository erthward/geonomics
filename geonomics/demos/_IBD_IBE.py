#!/usr/bin/python
# IBD_IBE_test.py

# flake8: noqa


########
# set-up
########

# geonomics imports
from geonomics.utils.viz import _check_display

# other imports
from copy import deepcopy
# from itertools import chain
import numpy as np
from sklearn.decomposition import PCA
import os
import matplotlib as mpl
_check_display()
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from copy import copy, deepcopy
import time
import statsmodels.api as sm




###############
# function defs
###############

def map_genetic_PCA(mod):
    """run genetic PCA and mapping individuals colored in RGB
    space by the first 3 PCs"""
    species = mod.comm[0]
    land = mod.land
    # get array of resulting genomic data (i.e. 'speciome'),
    # genotypes meaned by individual
    speciome = np.mean(np.stack([i.g for i in species.values()]), axis=2)
    # run PCA on speciome
    pca = PCA(n_components=3)
    PCs = pca.fit_transform(speciome)
    # normalize the PC results
    norm_PCs = (PCs - np.min(PCs,
                             axis=0)) / (np.max(PCs,
                                                axis=0) - np.min(PCs,
                                                                 axis=0))
    # use first 3 PCs to get normalized values for R, G, & B colors
    PC_colors = norm_PCs * 255
    # scatter all individuals on top of landscape, colored by the
    # RBG colors developed from the first 3 geonmic PCs
    xs = mod.comm[0]._get_x()
    ys = mod.comm[0]._get_y()
    # get environmental raster, with barrier masked out
    masked_env = deepcopy(mod.land[0].rast)
    masked_env[mod.land[1].rast == 0] = np.nan
    # create light colormap for plotting landscape
    # bot = plt.cm.get_cmap('Blues', 256)(np.linspace(0.4, 0.45, 2))[0]
    # top = plt.cm.get_cmap('Reds', 256)(np.linspace(0.4, 0.45, 2))[0]
    cmap = plt.cm.coolwarm
    copy(cmap).set_bad(color='#8C8C8C')
    # plot landscape
    # plt.imshow(masked_env, cmap=cmap, alpha=0.8)
    mod.plot(0, 0, color=PC_colors/255.0, mask_rast=masked_env,
             size=mark_size, edge_color='black')
    #plt.pcolormesh(land._x_cell_bds, land._y_cell_bds, masked_env, cmap=cmap)
    # scatter plot of individuals, colored by composite PC score
    #plt.scatter(xs, ys, c=PC_colors/255.0, s=mark_size, edgecolors='black')
    # fix x and y limits
    #[f([dim - 0.5 for dim in (0, mod.land.dim[0])]) for f in (plt.xlim,
    #                                                          plt.ylim)]
    # get rid of x and y ticks
    #[f([]) for f in (plt.xticks, plt.yticks)]


def plot_genetic_PCA(mod):
    """run a genetic PCA and plot individuals on first 2 PCs, colored by which
    side of the barrier they are on"""
    from copy import deepcopy
    from sklearn.decomposition import PCA
    figsize = 6
    species = mod.comm[0]
    # get array of resulting genomic data (i.e.
    # 'speciome'), genotypes meaned by individual
    speciome = np.mean(np.stack([i.g for i in species.values()]), axis=2)
    # run PCA on speciome
    pca = PCA(n_components=2)
    PCs = pca.fit_transform(speciome)
    # normalize the PC results
    norm_PCs = (PCs - np.min(PCs,
                    axis=0)) / (np.max(PCs, axis=0) - np.min(PCs, axis=0))
    # assign a value to each species, 0 or 1, indicating whether
    # they're located on the left or right vertical half of the landscape
    ind_colors = [0 if (
                ind.x < mod.land.dim[0]/2) else 1 for ind in species.values()]
    # plot individuals on PCs 1 and 2, colored by their landscape half
    fig = plt.figure(figsize=(figsize, figsize), dpi= 80,
                     facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    PC_dict = dict(zip(range(len(ind_colors)), ind_colors))
    left = [k for k,v in PC_dict.items() if v == 0]
    right = [k for k,v in PC_dict.items() if v == 1]
    patches = []
    patches.append(ax1.scatter(norm_PCs[left, 0], norm_PCs[left, 1],
                               color = '#00ffff'))
    patches.append(ax1.scatter(norm_PCs[right, 0], norm_PCs[right, 1],
                               color = '#ff00ff'))
    ax1.set_xlabel('genetic PC 1')
    ax1.set_ylabel('genetic PC 2')
    ax1.legend(patches, ['left of barrier', 'right of barrier'])


# calculate euclidean distance from two n-length vectors
def calc_euc(x, y):
    euc = np.sqrt(sum([(n-m)**2 for n, m in zip(x, y)]))
    return euc


# calculate lower-triangular of PCA-based Euclidean genetic distances between
# all individuals, using a 'speciome' (2d array of all individs' genomes)
def calc_dists(species, dist_type='gen', env_lyrs=None, return_flat=True,
               allele_freq_diff=False, biallelic=False):
    # calculate genetic distance as the euclidean distance between individuals
    # in genetic PC space
    if dist_type == 'gen':
        speciome = species._get_genotypes(biallelic=biallelic)
        if not allele_freq_diff:
            pca = PCA()
            vals = pca.fit_transform(speciome)
        else:
            vals = speciome
    # calculate geographic distance as the linear euclidean distance between
    # individuals
    elif dist_type == 'geo':
        vals = np.stack([np.array((i.x, i.y)) for i in species.values()])
    # calculate environmental distance as the euclidean distance between
    # individuals' environments, for all environmental layers specified by
    # the 'env_lyrs' argument
    elif dist_type == 'env':
        vals = np.stack([np.array(i.e)[env_lyrs] for i in species.values()])
    #calculate the phenotypic distances
    elif dist_type == 'phn':
        vals = np.stack([np.array(i.z)[env_lyrs] for i in species.values()])
    # print(vals)
    # print(vals.shape)
    n_ind = vals.shape[0]
    dist_mat = np.ones([n_ind] * 2) * np.nan
    # dist_vals = [[calc_euc(i, j) for j in vals[n:,
    #                                        :]] for n, i in enumerate(vals)]
    for i in range(n_ind):
        for j in range(0, i+1):
            if dist_type == 'gen' and allele_freq_diff:
                dist_mat[i, j] = np.sum(np.abs(
                                        vals[:, i] - vals[:, j])) / vals.shape[1]
            else:
                dist_mat[i, j] = calc_euc(vals[i, :], vals[j, :])
    # check that all diagonal values are 0
    assert np.all(np.diag(dist_mat) == 0), "Not all diagonal values are 0!"
    if dist_type == 'gen' and allele_freq_diff:
        assert (np.all(0 <= dist_mat[np.invert(np.isnan(dist_mat))])
                and np.all(dist_mat[np.invert(np.isnan(dist_mat))] <= 1)), (
            "You are calculating genetic distances using allele-frequency "
            "differences, which should be normalized to 0 <= diff <= 1, but "
            "you have non-NaN values outside that range!")

    if return_flat:
        # flatten the lower triangle, if return 1-d of values
        dists = dist_mat[np.tril_indices(dist_mat.shape[0], -1)]
        assert dists.size == (n_ind**2 - n_ind)/2, ("Length not equal "
                                                    "to n(n-1)/2!")
    else:
        # make it a symmetric dist matrix, if returning the matrix
        dist_mat[np.isnan(dist_mat)] = 0
        dists = dist_mat + dist_mat.T
        assert dists.size == n_ind**2, "Size not equal to n*n!"
    # dist_vals = [item[1:] for item in dist_vals]
    # dists = [*chain.from_iterable(dist_vals)]
    # assert that the length is correct
    return dists


# define a function to calculate the mean difference between phenotype
# and environment for a species
def calc_mean_z_e_diff(spp, trait_num=0):
    zs = spp._get_z().ravel()
    es = spp._get_e(lyr_num=spp.gen_arch.traits[trait_num].lyr_num)
    mean_diff = np.mean(np.abs(zs - es))
    return mean_diff


# define a function to calculate individuals' mean fitness
def calc_mean_fitness(spp):
    return np.mean(spp._calc_fitness())


def track_horiz_crossing(mod, zone_edges, tracker, count):
    """track all horizontal crossing of some vertical zone on the landscape,
    delineated by x-axis bounds in zone_edges argument, so that the model can
    be run with and without the barrier mask, to quantify the extent to which
    the barrier increases landscape resistance
    """
    for idx, ind in mod.comm[0].items():
        # individ is assigned value 0 if on left side of the crossing zone...
        if ind.x <= zone_edges[0]:
            curr_side = 0
        #...or 1 if on right side...
        elif ind.x >= zone_edges[1]:
            curr_side = 1
        # or else nan if they are inside the zone
        else:
            curr_side = np.nan
        # if individ is in the tracker then they were alive last time step,
        # so grab and compare their previous value to this time step's value
        if idx in tracker:
            diff = curr_side - tracker[idx]
            # ignore nans, becuase they are currently inside the zone, so
            # whether or not they are crossing cannot be decided
            if np.isnan(diff):
                pass
            # increment the crossing counter, if a crossing has occurred
            elif diff != 0:
                count +=1
        # now update the individ's value in the tracker (unless they have a
        # nan, because they're currently inside the zone, so maintain the
        # value from their last time outside the zone until they exit
        # the zone and it can be determined whether or not they crossed
        if not np.isnan(curr_side):
            tracker[idx] = curr_side
    return count


# function for setting a precise number of axis ticks, rounded to a fixed
# number of digits
def fix_yax_n_ticks_digits(ax, vals, n_ticks, n_digits):
    minval = min(vals)
    maxval = max(vals)
    tick_vals = np.linspace(minval, maxval, n_ticks)
    tick_labs = [('%%0.%if' % n_digits) % v for v in tick_vals]
    ax.set_yticks(tick_vals)
    ax.set_yticklabels(tick_labs)


# set some plotting params
img_dir = ('/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/'
           'img/final/')
ax_fontdict = {'fontsize': 25,
               'name': 'Bitstream Vera Sans'}
ttl_fontdict = {'fontsize': 25,
                'name': 'Bitstream Vera Sans'}
mark_size = 60


def _make_params():

    # IBD_IBE_params.py

    # This is a parameters file generated by Geonomics
    # (by the gnx.make_parameters_file() function).


                       ##  ::::::          :::    :: ::::::::::##
                 ##:::::    ::::   :::      ::    :: :: ::::::::::: :##
              ##::::::::     ::            ::   ::::::::::::::::::::::::##
            ##:::::::::                      :::::::::: :::::: ::::::::  :##
          ## : ::::  ::                    ::::  : ::    :::::::: : ::  :   ##
         ##GGGGG  EEEEE OOOOO   NN   NN   OOOOO   MM   MM IIIIII  CCCCC SSSSS##
        ##GG     EE    OO   OO  NNN  NN  OO   OO  MM   MM   II   CC     SS    ##
        ##GG     EE   OO     OO NN N NN OO     OO MMM MMM   II   CC     SSSSSS##
        ##GG GGG EEEE OO     OO NN  NNN OO     OO MM M MM   II   CC         SS##
        ##GG   G EE    OO   OO  NN   NN  OO   OO  MM   MM   II   CC        SSS##
         ##GGGGG  EEEEE OOOOO   NN   NN   OOOOO   MM   MM IIIIII  CCCCC SSSSS##
          ##    :::::::::               :::::::::: ::              ::  :   :##
            ##:   :::::                    :::::: :::             ::::::: ##
              ##   :::                      :::::  ::              :::::##
                 ## ::                      ::::                     ##
                       ##                                      ##
                          ## :: ::    :::            ##


    env_left_half = np.hstack([np.atleast_2d(np.linspace(
                    0, 1, 40) + np.random.normal(
                    0, 0.05, 40)).T for _ in range(20)])
    env_right_half = np.hstack([np.atleast_2d(np.linspace(
                    0, 1, 40) + np.random.normal(
                    0, 0.05, 40)).T for _ in range(20)])
    env_right_half = np.flipud(env_right_half)
    env = np.hstack((env_left_half, env_right_half))
    env = np.clip(env, 0, 1)
    barrier = np.ones((40, 40))
    barrier[:, 18:22] = 0


    params = {
    ###############################################################################

    ###################
    #### LANDSCAPE ####
    ###################
        'landscape': {

        ##############
        #### main ####
        ##############
            'main': {
                #y,x (a.k.a. i,j) dimensions of the Landscape
                'dim':                      (40,40),
                #resolution of the Landscape
                'res':                      (1,1),
                #upper-left corner of the Landscape
                'ulc':                      (0,0),
                #projection of the Landscape
                'prj':                      None,
                }, # <END> 'main'

        ################
        #### layers ####
        ################
            'layers': {

                #layer name (LAYER NAMES MUST BE UNIQUE!)
                'env': {

            #######################################
            #### layer num. 0: init parameters ####
            #######################################

                    #initiating parameters for this layer
                    'init': {
                        #parameters for a 'defined'-type Layer
                        'defined': {
                            'rast': env,
                            'pts': None,
                            'vals': None,
                            'interp_method': None
                            }, # <END> 'defined'

                        }, # <END> 'init'

                    }, # <END> layer num. 0

                #layer name (LAYER NAMES MUST BE UNIQUE!)
                'barrier': {

            #######################################
            #### layer num. 1: init parameters ####
            #######################################

                    #initiating parameters for this layer
                    'init': {
                        #parameters for a 'defined'-type Layer
                        'defined': {
                            'rast': barrier,
                            'pts': None,
                            'vals': None,
                            'interp_method': None
                            }, # <END> 'defined'

                        }, # <END> 'init'

                    }, # <END> layer num. 0



        #### NOTE: Individual Layers' sections can be copy-and-pasted (and
        #### assigned distinct keys and names), to create additional Layers.


                } # <END> 'layers'

            }, # <END> 'landscape'


    ###############################################################################

    ###################
    #### COMMUNITY ####
    ###################
        'comm': {

            'species': {

                #species name (SPECIES NAMES MUST BE UNIQUE!)
                'spp_0': {

                #####################################
                #### spp num. 0: init parameters ####
                #####################################

                    'init': {
                        #starting number of individs
                        'N':                1000,
                        #carrying-capacity Layer name
                        'K_layer':          'barrier',
                        #multiplicative factor for carrying-capacity layer
                        'K_factor':         1.5
                        }, # <END> 'init'

                #######################################
                #### spp num. 0: mating parameters ####
                #######################################

                    'mating'    : {
                        #age(s) at sexual maturity (if tuple, female first)
                        'repro_age':                0,
                        #whether to assign sexes
                        'sex':                      False,
                        #ratio of males to females
                        'sex_ratio':                1/1,
                        #intrinsic growth rate
                        'R':                        0.5,
                        #intrinsic birth rate (MUST BE 0<=b<=1)
                        'b':                        0.5,
                        #expectation of distr of n offspring per mating pair
                        'n_births_distr_lambda':    1,
                        #whether n births should be fixed at n_births_dist_lambda
                        'n_births_fixed':           True,
                        #radius of mate-search area
                        'mating_radius':            2,
                        'inverse_dist_mating':      False,
                        'choose_nearest_mate':      False,
                        }, # <END> 'mating'

                ##########################################
                #### spp num. 0: mortality parameters ####
                ##########################################

                    'mortality'     : {
                        #maximum age
                        'max_age':                      None,
                        #min P(death) (MUST BE 0<=d_min<=1)
                        'd_min':                        0,
                        #max P(death) (MUST BE 0<=d_max<=1)
                        'd_max':                        1,
                        #width of window used to estimate local pop density
                        'density_grid_window_width':    None,
                        }, # <END> 'mortality'

                #########################################
                #### spp num. 0: movement parameters ####
                #########################################

                    'movement': {
                        #whether or not the species is mobile
                        'move':                     True,
                        #mode of distr of movement direction
                        'direction_distr_mu':       1,
                        #concentration of distr of movement direction
                        'direction_distr_kappa':    0,
                        #mean of distr of movement distance
                        'movement_distance_distr_param1':        0.5,
                        #variance of distr of movement distance
                        'movement_distance_distr_param2':     0.5,
                        'movement_distance_distr':            'wald',
                        #first param of distr of dispersal distance
                        'dispersal_distance_distr_param1':       0.5,
                        #second param of distr of dispersal distance
                        'dispersal_distance_distr_param2':    0.5,
                        'dispersal_distance_distr':  'wald',
                        'move_surf'     : {
                            #move-surf Layer name
                            'layer': 'barrier',
                            #whether to use mixture distrs
                            'mixture': True,
                            #concentration of distrs
                            'vm_distr_kappa': 12,
                            #length of approximation vectors for distrs
                            'approx_len': 5000,
                            },  # <END> 'move_surf'

                        },    # <END> 'movement'


                #####################################################
                #### spp num. 0: genomic architecture parameters ####
                #####################################################

                    'gen_arch': {
                        #file defining custom genomic arch
                        'gen_arch_file':            None,
                        #num of loci
                        'L':                        100,
                        #num of chromosomes
                        'l_c':                      [100],
                        #whether starting allele frequencies should be fixed at 0.5
                        'start_p_fixed':            0.5,
                        #genome-wide per-base neutral mut rate (0 to disable)
                        'mu_neut':                  0,
                        #genome-wide per-base deleterious mut rate (0 to disable)
                        'mu_delet':                 0,
                        #shape of distr of deleterious effect sizes
                        'delet_alpha_distr_shape':  0.2,
                        #scale of distr of deleterious effect sizes
                        'delet_alpha_distr_scale':  0.2,
                        #alpha of distr of recomb rates
                        'r_distr_alpha':            None,
                        #beta of distr of recomb rates
                        'r_distr_beta':             None,
                        #whether loci should be dominant (for allele '1')
                        'dom':                      False,
                        #whether to allow pleiotropy
                        'pleiotropy':               False,
                        #custom fn for drawing recomb rates
                        'recomb_rate_custom_fn':    None,
                        #number of recomb paths to hold in memory
                        'n_recomb_paths_mem':       int(1e4),
                        #total number of recomb paths to simulate
                        'n_recomb_paths_tot':       int(1e5),
                        'n_recomb_sims':            10000,
                        'start_neut_zero':          False,
                        'allow_ad_hoc_recomb':      False,
                        #whether to save mutation logs
                        'mut_log':                  False,

                        'traits': {

                            ###########################
                            ####trait 0 parameters ####
                            ###########################
                            #trait name (TRAIT NAMES MUST BE UNIQUE!)
                            'trait_0': {
                                #trait-selection Layer name
                                'layer':                'env',
                                #polygenic selection coefficient
                                'phi':                  0.05,
                                #number of loci underlying trait
                                'n_loci':               10,
                                #mutation rate at loci underlying trait
                                'mu':                   0,
                                #mean of distr of effect sizes
                                'alpha_distr_mu':      0.1,
                                #variance of distr of effect size
                                'alpha_distr_sigma':    0,
                                #max alpha value
                                'max_alpha_mag':        None,
                                #curvature of fitness function
                                'gamma':                1,
                                #whether the trait is universally advantageous
                                'univ_adv':             False
                                }, # <END> trait 0


        #### NOTE: Individual Traits' sections can be copy-and-pasted (and
        #### assigned distinct keys and names), to create additional Traits.


                            }, # <END> 'traits'

                        }, # <END> 'gen_arch'


                    }, # <END> spp num. 0



        #### NOTE: individual Species' sections can be copy-and-pasted (and
        #### assigned distinct keys and names), to create additional Species.


                }, # <END> 'species'

            }, # <END> 'comm'


    ###############################################################################

    ###############
    #### MODEL ####
    ###############
        'model': {
            #total Model runtime (in timesteps)
            'T':            10,
            #min burn-in runtime (in timesteps)
            'burn_T':       30,
            #seed number
            'num':          None,
            #time step interval for simplification of tskit tables
            'tskit_simp_interval':      100,

            ###############################
            #### iterations parameters ####
            ###############################
            'its': {
                #num iterations
                'n_its':            1,
                #whether to randomize Landscape each iteration
                'rand_landscape':   False,
                #whether to randomize Community each iteration
                'rand_comm':        False,
                #whether to burn in each iteration
                'repeat_burn':      False,
                }, # <END> 'iterations'



            } # <END> 'model'

        } # <END> params

    return params


def _run(params, save_figs=False, time_it=False,
         make_3d_plots=False, use_barrier=True,
         run_mantel=False):
    """
    use_barrier:
        sets flag to indicate whether or not to use the barrier
        (this will typically be used, but can be turned off to run the model
        without a barrier and thus determine, comparatively, how much the
        barrier increases landscape resistance)
    """

    ###############
    # set up figure
    ###############
    nrows = 13
    ncols = 4
    ratio = ncols/nrows
    vert_size = 11 # length of figure's vertical dimension
    #fig = plt.figure(figsize=(vert_size, vert_size*ratio))
    plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.96, wspace=0.07,
                        hspace=0.16)
    gs = gridspec.GridSpec(nrows, ncols)
                          #height_ratios=[1]*8 + [1.5]*6,
                          #width_ratios=[1,1,1,0.7,0.7,1,1,1])

    #fig bounds
    gen_b4_top = 0
    gen_b4_bot = 3
    gen_b4_L = 0
    gen_b4_R = 2

    gen_af_top = 0
    gen_af_bot = 3
    gen_af_L = 2
    gen_af_R = 4

    ze_fit_top = 6
    ze_fit_bot = 10
    ze_fit_L = 1
    ze_fit_R = 3

    phn_b4_top = 3
    phn_b4_bot = 6
    phn_b4_L = 0
    phn_b4_R = 2

    phn_af_top = 3
    phn_af_bot = 6
    phn_af_L = 2
    phn_af_R = 4

    n1_3d_top = 10
    n1_3d_bot = 13
    n1_3d_L = 0
    n1_3d_R = 2

    n3_3d_top = 10
    n3_3d_bot = 13
    n3_3d_L = 2
    n3_3d_R = 4

    # BEFORE-SIM AXES
    gen_b4_ax = plt.subplots(1,1)[1]
    #gen_b4_ax = fig.add_subplot(gs[gen_b4_top:gen_b4_bot,
    #                               gen_b4_L:gen_b4_R], aspect='equal')
    gen_b4_ax.set_ylabel('genotype', fontdict=ax_fontdict)
    gen_b4_ax.set_title('before simulation', fontdict=ttl_fontdict)
    phn_b4_ax = plt.subplots(1,1)[1]
    #phn_b4_ax = fig.add_subplot(gs[phn_b4_top:phn_b4_bot,
    #                               phn_b4_L:phn_b4_R], aspect='equal')
    phn_b4_ax.set_ylabel('phenotype', fontdict=ax_fontdict)

    # AFTER-SIM AXES
    gen_af_ax = plt.subplots(1,1)[1]
    #gen_af_ax = fig.add_subplot(gs[gen_af_top:gen_af_bot,
    #                               gen_af_L:gen_af_R], aspect='equal')
    gen_af_ax.set_title('after simulation', fontdict=ttl_fontdict)
    #gen_af_ax.set_ylabel('genotype', fontdict=ax_fontdict)
    phn_af_ax = plt.subplots(1,1)[1]
    #phn_af_ax = fig.add_subplot(gs[phn_af_top:phn_af_bot,
    #                               phn_af_L:phn_af_R], aspect='equal')
    #phn_af_ax.set_ylabel('phenotype', fontdict=ax_fontdict)

    #--------
    # 3D AXES
    #--------
    # num 1
    n1_3d_ax = plt.figure().add_subplot(111, projection='3d')
    #n1_3d_ax = fig.add_subplot(gs[n1_3d_top:n1_3d_bot,
                                  #n1_3d_L:n1_3d_R], projection='3d')
    n1_3d_ax.view_init(elev=3, azim=83)
    n1_3d_ax.set_xlabel('$\longleftarrow$ geo. dist.', size=25, labelpad=-13)
    #n1_3d_ax.set_ylabel(' ' * 35 + '$\longleftarrow$ Env. Dist.', size=9,
    #                    labelpad=20)
    n1_3d_ax.zaxis.set_rotate_label(False)
    n1_3d_ax.set_zlabel('gen. dist. $\longrightarrow$', size=25, labelpad=-13,
                        rotation=90)
    n1_3d_ax.set_xticklabels([])
    n1_3d_ax.set_yticklabels([])
    n1_3d_ax.set_zticklabels([])
    n1_3d_ax.set_title("", pad=-260)
    # num 2
    ze_fit_ax_L = plt.subplots(1,1)[1]
    #ze_fit_ax_L = fig.add_subplot(gs[ze_fit_top:ze_fit_bot,
    #                                 ze_fit_L:ze_fit_R])
    ze_fit_ax_R = ze_fit_ax_L.twinx()
    #n2_3d_ax = fig.add_subplot(gs[2, 3:6], projection='3d')
    #n2_3d_ax.view_init(elev=25, azim=225)
    #n2_3d_ax.set_xlabel('Geo. Dist. $\longrightarrow$', size=9, labelpad=-13)
    #n2_3d_ax.set_ylabel('$\longleftarrow$ Env. Dist.', size=9, labelpad=-13)
    #n2_3d_ax.zaxis.set_rotate_label(False)
    #n2_3d_ax.set_zlabel('Gen. Dist. $\longrightarrow$', size=9, labelpad=-13,
    #                    rotation=90)
    #n2_3d_ax.set_xticklabels([])
    #n2_3d_ax.set_yticklabels([])
    #n2_3d_ax.set_zticklabels([])

    # num 3
    n3_3d_ax = plt.figure().add_subplot(111, projection='3d')
    #n3_3d_ax = fig.add_subplot(gs[n3_3d_top:n3_3d_bot,
    #                              n3_3d_L:n3_3d_R], projection='3d')
    n3_3d_ax.view_init(elev=3, azim=7)
    #n3_3d_ax.set_xlabel('$\longleftarrow$ Geo. Dist.' + ' ' * 25, size=9,
    #                    labelpad=10)
    n3_3d_ax.set_ylabel('env. dist. $\longrightarrow$', size=25, labelpad=-13)
    n3_3d_ax.zaxis.set_rotate_label(False)
    #n3_3d_ax.set_zlabel('gen. dist. $\longrightarrow$', size=25, labelpad=-13,
    #                    rotation=90)
    n3_3d_ax.set_xticklabels([])
    n3_3d_ax.set_yticklabels([])
    n3_3d_ax.set_zticklabels([])



    #########################
    # set trackers and timers
    #########################

    # create objects for crossing of barrier, and of equal-width vertical areas
    # within each of the two sides for comparison
    cross_count = 0
    cross_tracker = {}

    # start timer
    if time_it:
        start = time.time()


    #####################
    # prep and make model
    #####################

    # set model name (since the params are'nt being read in from separate file)
    params.model['name'] = 'IBD_IBE_demo'

    # get barrier-zone edges, for tracking crossings
    barr_rast = params['landscape']['layers']['barrier']['init']['defined']['rast']
    zone_edges = np.where(barr_rast[0,:] == 0)[0]
    zone_edges = (zone_edges.min(), zone_edges.max() + 1)

    # then get rid of barrier, if not to be used for this run
    if not use_barrier:
        params['landscape']['layers']['barrier']['init']['defined']['rast'] = np.ones(
                                            params['landscape']['main']['dim'])

    from .. import make_model
    mod = make_model(params)




    ################################
    # run model, plotting as it goes
    ################################

    # define number of timesteps
    T = 1000

    # burn model in
    mod.walk(20000, 'burn')

    # plot genetic PCA before evolution begins
    plt.sca(gen_b4_ax)
    map_genetic_PCA(mod)

    # plot phenotypes before evolution begins
    mask = np.ma.masked_where(mod.land[1].rast == 0, mod.land[1].rast)
    plt.sca(phn_b4_ax)
    mod.plot_phenotype(0, 0, mask_rast=mask, size=mark_size, cbar=False)
    [f((-0.5, 39.5)) for f in [phn_b4_ax.set_xlim, phn_b4_ax.set_ylim]]

    # create data structure to save z-e diff values
    if not time_it:
        mean_z_e_diffs = []
        mean_fits = []

    # run model for T timesteps
    for t in range(T):
        if not time_it:
            mean_z_e_diffs.append(calc_mean_z_e_diff(mod.comm[0]))
            mean_fits.append(calc_mean_fitness(mod.comm[0]))
            cross_count = track_horiz_crossing(mod, zone_edges, cross_tracker,
                                               cross_count)
        mod.walk(1)
    if not time_it:
        mean_z_e_diffs.append(calc_mean_z_e_diff(mod.comm[0]))
        mean_fits.append(calc_mean_fitness(mod.comm[0]))
        cross_count = track_horiz_crossing(mod, zone_edges, cross_tracker,
                                           cross_count)

    # finish timer
    if time_it:
        stop = time.time()
        tot_time = stop - start


    # plot genetic PCA after 1/4T timesteps
    plt.sca(gen_af_ax)
    map_genetic_PCA(mod)

    # plot the individuals' phenotypes
    plt.sca(phn_af_ax)
    mod.plot_phenotype(0, 0, mask_rast=mask, size=mark_size, cbar=False)
    [f((-0.5, 39.5)) for f in [phn_af_ax.set_xlim, phn_af_ax.set_ylim]]



    #########################
    # create 3d IBD, IBE plot
    #########################

    spp_subset = deepcopy(mod.comm[0])
    rand_inds = spp_subset._get_random_individuals(250)
    all_inds = [*spp_subset]
    for ind in all_inds:
        if ind not in rand_inds:
            spp_subset.pop(ind)
    gen_dists = calc_dists(spp_subset, biallelic=False, allele_freq_diff=False)
    # output the genetic distance matrix to CSV, if Mantel test required
    if run_mantel:
        gen_dist_mat = calc_dists(spp_subset, biallelic=False,
                                  allele_freq_diff=False,return_flat=False)
        np.savetxt("gen.csv", gen_dist_mat, delimiter=',')
    scaled_gen_dists = gen_dists/gen_dists.max()
    assert (np.all(scaled_gen_dists >= 0)
            and np.all(scaled_gen_dists <= 1)), ('Scaled genetic dist is outside '
                                                 '0 and 1!')
    geo_dists = calc_dists(spp_subset, 'geo')
    env_dists = calc_dists(spp_subset, 'env', [0])
    pheno_dists = calc_dists(spp_subset, 'phn', [0])

    # output the geo and env distance matrices to CSV, if Mantel test required,
    # then run the Mantel test
    if run_mantel:
        geo_dist_mat = calc_dists(spp_subset, 'geo', return_flat=False)
        env_dist_mat = calc_dists(spp_subset, 'env', [0], return_flat=False)
        np.savetxt("geo.csv", geo_dist_mat, delimiter=',')
        np.savetxt("env.csv", env_dist_mat, delimiter=',')
        # get the path to the mantel script, then run it
        gnx_path = '/' + os.path.join(*os.path.split(__file__)[0].split('/')[:-1])
        mantel_script_dir = os.path.join(gnx_path, 'data/IBD_IBE_demo')
        mantel_script = os.path.join(mantel_script_dir,
                                     os.listdir(mantel_script_dir)[0])
        print('\n\n')
        os.system('Rscript %s' % mantel_script)
        print('\n\n')

    # run multiple linear regression of gen on geo and env dists
    mlr_est = sm.Logit(endog=np.array(scaled_gen_dists.T),
                       exog=np.vstack((geo_dists, env_dists)).T).fit()

    # create logistic regression (to use in wireframe predicted surface)
    x_preds = np.arange(0, 50.1, 0.1)
    y_preds = np.linspace(0, 1, len(x_preds))
    z_preds = mlr_est.predict(np.vstack((x_preds, y_preds)).T)
    p_val_lt = mlr_est.pvalues[0] < 0.001
    assert p_val_lt, 'p-value not less than 0.001!'
    print(mlr_est.summary2())

    # create 3d-plot data
    y_vals = np.arange(0, 1.01, 0.01)
    ys = np.hstack([list(y_vals) for _ in range(len(y_vals))])
    xs = np.hstack([[n] * len(y_vals) for n in np.linspace(0, 50, len(y_vals))])
    zs = mlr_est.predict(np.vstack((xs, ys)).T)
    xs = xs.reshape([len(y_vals)] * 2)
    ys = ys.reshape([len(y_vals)] * 2)
    zs = zs.reshape([len(y_vals)] * 2)
    col3d = pheno_dists / pheno_dists.max()


    # plot on 3d axes
    for ax in [n1_3d_ax, n3_3d_ax]:#, n2_3d_ax]:
       ax.scatter(geo_dists, env_dists, scaled_gen_dists,
                  alpha=0.7, edgecolor='black', c=col3d, cmap='plasma')

    ##########################
    # create plot of z-e diffs
    ##########################
    if not time_it:
        L_color = '#096075'
        R_color = '#bf2659'
        ze_fit_ax_L.set_xlabel('time (steps)', fontdict=ax_fontdict)
        ze_fit_ax_L.set_ylabel('mean($|z-e|$)', fontdict=ax_fontdict,
                               color=L_color)
        ze_fit_ax_L.tick_params(axis='y', labelcolor=L_color, labelrotation=45)
        ze_fit_ax_L.plot(range(len(mean_z_e_diffs)), mean_z_e_diffs, color=L_color)
        fix_yax_n_ticks_digits(ze_fit_ax_L, mean_z_e_diffs, 5, 2)
        ze_fit_ax_R.set_ylabel('mean fitness', color=R_color,
                               fontdict=ax_fontdict)
        ze_fit_ax_R.tick_params(axis='y', labelcolor=R_color, labelrotation=-45)
        ze_fit_ax_R.plot(range(len(mean_fits)), mean_fits, color=R_color)
        fix_yax_n_ticks_digits(ze_fit_ax_R, mean_fits, 5, 3)
        #ax.set_xlabel('time')
        #ax.set_ylabel(('mean difference between individuals\' phenotypes and '
        #               'environmental values'))
        #z_e_fig.show()
        #if save_figs:
        #    fig.savefig('IBD_IBE_z-e_plot.png', format='png', dpi=1000)

        #fig.tight_layout()
        #fig.show()
        if save_figs:
                fig.savefig('IBD_IBE.png', format='png', dpi=1000)


    #########################
    # create 3d animated plot
    #########################
    if make_3d_plots:
        plt.rc('animation', html='html5')
        fig3d = plt.figure()
        ax3d = fig3d.add_subplot(111, projection='3d')
        ax3d.set_xlabel('geo', size=15)
        ax3d.set_ylabel('env', size=15)
        ax3d.set_zlabel('gen', size=15)
        # initialization function, which plots the background of each frame
        def init():
            ax3d.plot_wireframe(xs, ys, zs, color='lightgray', ccount=10, rcount=10,
                                linestyles='dashed', alpha=0.95, linewidths=0.5)
            scat = ax3d.scatter(geo_dists, env_dists, scaled_gen_dists,
                                 alpha=0.5, c=col3d, cmap='plasma')
            return fig3d,

        # animation function, which is called sequentially
        def animate(i):
            ax3d.view_init(elev=5., azim=i)
            return fig3d,

        # use the init and animate functions to create an animation object
        anim = animation.FuncAnimation(fig3d, animate, init_func=init, frames=359,
                                       interval=5, blit=True)
        # write to file
        if save_figs:
            try:
                anim.save('IBD_IBE_animation.gif', writer='imagemagick',
                          fps=60)
            except Exception as e:
                print(('\nCould not use Imagemagick to save the 3D plot\'s '
                       'GIF. The following error was thrown:\n\t%s') % e)

        fig3d.show()


    #####################################
    # print out the zone-crossing results
    #####################################
    cross_rate = cross_count / T
    print(("\nThe barrier-zone crossing rate was %0.6f individuals "
          "per time step.\n") % cross_rate)

    # print out time
    if time_it:
        print("\n\nModel ran in %0.2f seconds." % tot_time)

    return mod

