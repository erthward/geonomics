#!/usr/bin/python

#####################
#TODO
#1.  Add a "timing" argument, that I can set to false to run the model without
#plotting, to get an accurate assessment of run time

#2. Debug all the plotting stuff
#####################

# geonomics imports
from geonomics.utils.viz import _check_display

# other imports
import os
import numpy as np
import matplotlib as mpl
_check_display()
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d.axes3d import Axes3D
from nlmpy import nlmpy
from pykrige.ok import OrdinaryKriging
import time
import rasterio as rio
from scipy.stats import norm


gnx_dir = os.path.split(__file__)[0]
gnx_dir = os.path.join(*os.path.split(gnx_dir)[:-1])
DATA_PATH = os.path.join(gnx_dir, "data", "yosemite_demo")


def calc_neighborhood_mean_phenotype(mod, window_width=8,
                                     rel_bandwidth=0.5):
    kde = True
    if kde:
        print('Now calculating phenotype KDE. Please wait...\n\n')
        xs = mod.comm[0]._get_x()
        ys = mod.comm[0]._get_y()
        zs = mod.comm[0]._get_z()[:,0]
        OK = OrdinaryKriging(xs, ys, zs, variogram_model='spherical')
        gridx = np.arange(0.5, mod.land.dim[0] + 0.5, 1)
        gridy = np.arange(0.5, mod.land.dim[1] + 0.5, 1)
        z, ss = OK.execute('grid', gridx, gridy)
        return z
    if not kde:
        bandwidth = window_width * rel_bandwidth
        # array to store each cell's mean phenotype value
        mean_z = np.ones(mod.land.dim)
        # calc half-window_width
        hww = int(window_width / 2)
        # create the kernel to use for weighting individuals' phenotypes
        pd = norm(0, bandwidth)
        # calculate the max prob dens (to normalize weights in the weights lists
        # below)
        max_weight = pd.pdf(0)
        print('max weight is', max_weight)
        # loop over cells
        for i in range(mod.land.dim[1]):
            for j in range(mod.land.dim[0]):
                # get window around each cell
                i_min = max(i + 0.5 - hww, 0)
                i_max = min(i + 0.5 + hww, mod.land.dim[1])
                j_min = max(j + 0.5 - hww, 0)
                j_max = min(j + 0.5 + hww, mod.land.dim[0])
                print(i_min, i_max, j_min, j_max)
                # get all phenotypes in the window, Gaussian-weighted by their
                # distances from the window center
                zs_in_window = []
                dists_from_center = []
                effective_zs = []
                for ind in mod.comm[0].values():
                    if ((i_min <= ind.y <= i_max)
                        and (j_min <= ind.x <= j_max)):
                        print(ind.x, ind.y, 'are inside the zone!')
                        zs_in_window.append(ind.z[0])
                        dist = np.sqrt((i + 0.5 - ind.x)**2 + (j + 0.5 - ind.y)**2)
                        dists_from_center.append(dist)
                # average the window's phenotypes and add to mean_z
                # NOTE: if there are no individuals in the window then a NaN is
                # returned
                print('ziw')
                print(zs_in_window)
                print('dfc')
                print(dists_from_center)
                weights = [pd.pdf(d)/max_weight for d in dists_from_center]
                print('w')
                print(weights)
                effective_zs = [zs_in_window[n] * weights[n] for n in range(
                                                                    len(weights))]
                print('ez')
                print(effective_zs)
                print('\n\n\n------------------\n\n\n')
                mean_z[i, j] = np.mean(effective_zs)

        return mean_z


# plot 3 change rasters, 1 for each of the 3 landscape layers
def make_change_plt(mod):
    change_fig = plt.figure()
    change_fig_titles = {0: 'temperature',
                         2: 'precipitation',
                         1: 'habitat suitability',
                     }
    cmaps = ['YlOrRd', 'Spectral', 'PiYG']
    for i, d in enumerate(change_fig_titles.items()):
        ax = change_fig.add_subplot(1, 3, i + 1)
        lyr_num = d[0]
        title = d[1]
        ax.set_title(title, size=20)
        af = mod.land[lyr_num]._get_rast_in_native_units()
        b4 = mod.orig_land[lyr_num]._get_rast_in_native_units()
        diff = af - b4
        im = ax.imshow(diff, cmap = cmaps[i])
        units = '$\Delta$ ' + mod.land[lyr_num].units
        cb = change_fig.colorbar(im, orientation='horizontal')
        cb.set_label(label=units, size=13)
        cb.ax.tick_params(labelsize=9)
        (xtick_locs, xtick_labs, ytick_locs,
         ytick_labs) = mod.land[lyr_num]._get_coord_ticks()
        ax.set_xticks(xtick_locs)
        ax.set_xticklabels(np.round(xtick_labs, 1), rotation=45, size=7,
                           color='gray')
        ax.set_yticks(ytick_locs)
        ax.set_yticklabels(ytick_labs, size=7, color='gray')
        ax.set_xlabel('lon', color='gray', size=10)
        if i == 0:
            ax.set_ylabel('lat', color='gray', size=10)
    change_fig.show()


# drape a raster on top of a DEM
def drape_raster(mod, rast, DEM, cmap='rainbow'):
    # get the x and y coordinates for each cell
    xi = range(mod.land.dim[0])
    yi = range(mod.land.dim[1])
    # use the x and y coords to get a meshgrid
    X, Y = np.meshgrid(xi, yi)
    # create a 3d-axes plot
    drape_fig = plt.figure()
    ax = Axes3D(drape_fig)
    # create a mappable object, to use to make the colorbar
    mappable = plt.cm.ScalarMappable()
    mappable.set_array(rast)
    mappable.set_cmap(cmap)
    # plot the DEM as a surface, coloring its patches by the draped raster
    surf = ax.plot_surface(X, Y, DEM, facecolors=mappable.cmap(rast),
                           cmap=mappable.cmap)
    # create a colorbar
    cb = plt.colorbar(mappable)
    cb.set_label(label='index of phenotypic change', size=13)
    cb.ax.tick_params(labelsize=9)
    # set the x-, y-, and z-axis ticks and labels
    (xtlocs, xtlabs, ytlocs, ytlabs) = mod.land[0]._get_coord_ticks()
    ztlocs = ax.zaxis.get_ticklocs()
    ztlabs = ax.zaxis.get_ticklabels()
    ax.set_xticks(xtlocs)
    ax.set_xticklabels(np.round(xtlabs, 1), rotation=45, size=7, color='gray')
    ax.set_yticks(ytlocs)
    ax.set_yticklabels(np.round(ytlabs, 1), size=7, color='gray')
    ax.set_zticks(ztlocs)
    ax.set_zticklabels(ztlabs, size=7, color='gray')
    ax.set_xlabel('lon', color='gray', size=10)
    ax.set_ylabel('lat', color='gray', size=10)
    ax.set_zlabel('alt (m)', color='gray', size=10)
    drape_fig.show()
    return drape_fig


# create and save images that Imagemagick will stitch into a GIF
def save_gif_img(t):
    # set up second figure, for a GIF animation
    mid_col_width = 0.8
    fig2 = plt.figure(figsize=(6.75 + mid_col_width * 6.75, 5.4))
    ax2_gs = fig.add_gridspec(1, 3, width_ratios=[1, scnd_col_width, 1])
    # fig2, axs2 = plt.subplots(1, 2)
    # fig2.set_tight_layout(True)
    ax2_1 = fig2.add_subplot(ax2_gs[0, 0])
    ax2_2 = fig2.add_subplot(ax2_gs[0, 2])
    fig2.suptitle('Time step: %i' % (t+1), size=20)
    ax2_1.set_title('Temperature', size=12)
    ax2_2.set_title('Habitat suitability', size=12)
    ax2_1.set_xlabel('lon', size=8)
    ax2_1.set_ylabel('lat', size=8)
    ax2_2.set_xlabel('lon', size=8)
    ax2_2.set_ylabel('lat', size=8)
    ax2_1.set_aspect(1)
    ax2_2.set_aspect(1)
    # plot the population on both the temp and hab rasters
    land1 = ax2_1.pcolormesh(mod.land._x_cell_bds,
                             mod.land._y_cell_bds,
                             mod.land[0]._get_rast_in_native_units(),
                             cmap='coolwarm',
                             vmin=mod.land[0]._scale_min,
                             vmax=mod.land[0]._scale_max)
    xticks, xticklabs, yticks, yticklabs = mod.land[0]._get_coord_ticks()
    ax2_1.set_xticklabels(labels=xticklabs, size=6)
    ax2_1.set_yticklabels(labels=yticklabs, size=6)
    divider = make_axes_locatable(ax2_1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar_ticks = mod.land[0]._get_cbar_ticks_and_minmax_scaled_vals()
    cbar = plt.colorbar(land1, cax=cax, ticks=np.linspace(cbar_ticks[1],
                                                          cbar_ticks[2], 5))
    cbar.ax.set_yticklabels(cbar_ticks[0], size=6)
    cbar.set_label(mod.land[0].units, rotation=270, labelpad=5, y=0.5,
                   size=7)
    coords = mod.comm[0]._get_plot_coords()
    xs = coords[:, 0]
    ys = coords[:, 1]
    zs = mod.comm[0]._get_z()[:, 0]
    ax2_1.scatter(xs, ys, c=zs, s=1, cmap='coolwarm')
    land2 = ax2_2.pcolormesh(mod.land._x_cell_bds,
                             mod.land._y_cell_bds,
                             mod.land[1]._get_rast_in_native_units(),
                             cmap='BrBG_r',
                             vmin=mod.land[1]._scale_min,
                             vmax=mod.land[1]._scale_max)
    xticks, xticklabs, yticks, yticklabs = mod.land[1]._get_coord_ticks()
    ax2_2.set_xticklabels(labels=xticklabs, size=6)
    ax2_2.set_yticklabels(labels=yticklabs, size=6)
    divider = make_axes_locatable(ax2_2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar_ticks = mod.land[1]._get_cbar_ticks_and_minmax_scaled_vals()
    cbar = plt.colorbar(land2, cax=cax, ticks=np.linspace(cbar_ticks[1],
                                                          cbar_ticks[2], 5))
    cbar.ax.set_yticklabels(cbar_ticks[0], size=6)
    cbar.set_label(mod.land[1].units, rotation=270, labelpad=15, y=0.5,
                   size=10)
    ax2_2.scatter(xs, ys, c='black', s=0.5)

    # save the figure
    if save_figs and make_gif:
        fig2.savefig('gif_img_%i.jpg' % (t + 1))

    # manually close the figure
    plt.close(fig2)



def _make_params():

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

                # NOTE: The resolution of the rasters is 0.00833 dec. deg.,
                # which is equal to about 0.00833* 88070.4 = 730.984 m in E-W
                # and about 0.00833 * 111320 = 927.2956 m in N-S directions.
                # NOTE: Thus, a cell is about 730.984 * 927.2956 = 678000 m^2,
                # which is = 67.8 hectares.
                # NOTE: The difference in distance in E-W and N-S is problematic
                # for distance-based things!!!!

                #y,x (a.k.a. i,j) dimensions of the Landscape
                'dim':                      (157, 157),
                #resolution of the Landscape
                'res':                      (0, 0),
                #upper-left corner of the Landscape
                'ulc':                      (0, 0),
                #projection of the Landscape
                'prj':                      None,
                }, # <END> 'main'

        ################
        #### layers ####
        ################
            'layers': {

                #layer name (LAYER NAMES MUST BE UNIQUE!)
                'tmp': {

            #######################################
            #### layer num. 0: init parameters ####
            #######################################

                    #initiating parameters for this layer
                    'init': {

                        #parameters for a 'file'-type Layer
                        'file': {
                            #</path/to/file>.<ext>
                            'filepath': os.path.join(DATA_PATH,
                                                     'yosemite_lyrs',
                                                     'tmp_1980-2010.tif'),
                            #minimum value to use to rescale the Layer to [0,1]
                            'scale_min_val':                -1.37,
                            #maximum value to use to rescale the Layer to [0,1]
                            'scale_max_val':                19.75773,
                            #decimal precision to use for coord-units (ulc & res)
                            'coord_prec':                   8,
                            #units of this file's variable
                            'units':                    '$^∘C$',

                            }, # <END> 'file'

                        }, # <END> 'init'

                #########################################
                #### layer num. 0: change parameters ####
                #########################################

                    #landscape-change events for this Layer
                    'change': {

                        0: {
                            #array or file for final raster of event, or directory
                            #of files for each stepwise change in event
                            'change_rast': os.path.join(DATA_PATH,
                                                        'yosemite_lyrs',
                                                        'tmp'),

                            #starting timestep of event
                            'start_t':          509,
                            #ending timestep of event
                            'end_t':            594,
                            #number of stepwise changes in event
                            'n_steps':          18,
                            }, # <END> event 0

                        }, # <END> 'change'

                    }, # <END> layer num. 0


                #layer name (LAYER NAMES MUST BE UNIQUE!)
                'hab': {

            #######################################
            #### layer num. 1: init parameters ####
            #######################################

                    #initiating parameters for this layer
                    'init': {

                        #parameters for a 'file'-type Layer
                        'file': {
                            #</path/to/file>.<ext>
                            'filepath': os.path.join(DATA_PATH,
                                                     'yosemite_lyrs',
                                                     'sdm_1980-2010.tif'),
                            #minimum value to use to rescale the Layer to [0,1]
                            'scale_min_val':                0,
                            #maximum value to use to rescale the Layer to [0,1]
                            'scale_max_val':                1,
                            #decimal precision to use for coord-units (ulc & res)
                            'coord_prec':                   8,
                            #units of this file's variable
                            'units':                    'suitability',

                            }, # <END> 'file'

                        }, # <END> 'init'

                #########################################
                #### layer num. 1: change parameters ####
                #########################################

                    #landscape-change events for this Layer
                    'change': {

                        0: {
                            #array or file for final raster of event, or directory
                            #of files for each stepwise change in event
                            'change_rast': os.path.join(DATA_PATH,
                                                        'yosemite_lyrs',
                                                        'sdm'),

                            #starting timestep of event
                            'start_t':          509,
                            #ending timestep of event
                            'end_t':            594,
                            #number of stepwise changes in event
                            'n_steps':          18,
                            }, # <END> event 0

                        }, # <END> 'change'

                    }, # <END> layer num. 1


                #layer name (LAYER NAMES MUST BE UNIQUE!)
                'ppt': {

            #######################################
            #### layer num. 2: init parameters ####
            #######################################

                    #initiating parameters for this layer
                    'init': {

                        #parameters for a 'file'-type Layer
                        'file': {
                            #</path/to/file>.<ext>
                            'filepath': os.path.join(DATA_PATH,
                                                     'yosemite_lyrs',
                                                     'ppt_1980-2010.tif'),
                            #minimum value to use to rescale the Layer to [0,1]
                            'scale_min_val':                81.53713,
                            #maximum value to use to rescale the Layer to [0,1]
                            'scale_max_val':                2171.341,
                            #decimal precision to use for coord-units (ulc & res)
                            'coord_prec':                   8,
                            #units of this file's variable
                            'units':                    '$mm/yr$',

                            }, # <END> 'file'

                        }, # <END> 'init'

                #########################################
                #### layer num. 0: change parameters ####
                #########################################

                    #landscape-change events for this Layer
                    'change': {

                        0: {
                            #array or file for final raster of event, or directory
                            #of files for each stepwise change in event
                            'change_rast': os.path.join(DATA_PATH,
                                                        'yosemite_lyrs',
                                                        'ppt'),

                            #starting timestep of event
                            'start_t':          509,
                            #ending timestep of event
                            'end_t':            594,
                            #number of stepwise changes in event
                            'n_steps':          18,
                            }, # <END> event 0

                        }, # <END> 'change'

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
                'Sceloporus graciosus': {

                #####################################
                #### spp num. 0: init parameters ####
                #####################################

                    'init': {
                        #starting number of individs
                        'N':                5000,
                        #carrying-capacity Layer name
                        'K_layer':          'hab',
                        #multiplicative factor for carrying-capacity layer
                        # NOTE: each cell is ~53.43 hectares, and pop density
                        # estimates in Utah (closest estimate I could find to
                        # Yosemite region) were 208/hectare, with little variation
                        # (Tinkle 1973).
                        # That would suggest each cell should have a carrying
                        # capacity of ~ 53.43 * 208. However, in the Sierra they 
                        # and S. occidentalis segregate by habitat (Rose 1976),
                        # so if we assume that only about 10% of land area
                        # is covered by S. graciosus' preferred open,
                        # exposed habitat, then that suggests we should
                        # use a K_factor of 67.8 * 208 * 0.1 = 1111.344
                        'K_factor':         3,
                        }, # <END> 'init'

                #######################################
                #### spp num. 0: mating parameters ####
                #######################################

                    'mating'    : {
                        #age(s) at sexual maturity (if tuple, female first)
                        # NOTE: average reproductive age is 2 (i.e. in the third
                        # season) (Tinkle 1973; Tinkle et. al 1993)
                        'repro_age':                2,
                        #whether to assign sexes
                        'sex':                      True,
                        # NOTE: Some studies suggest a skewed sex ratio
                        # becuase of uneven male/female survival rate ratio
                        # (e.g. Tinkle 1973), but other, more recent studies do not
                        # (e.g. Tinkle et al. 1993), so just leaving sex_ratio at 1
                        #ratio of males to females
                        'sex_ratio':                1/1,
                        #whether P(birth) should be weighted by parental dist
                        'dist_weighted_birth':      False,
                        # NOTE: found no information about this, but 0.5 seems
                        # a reasonable default value
                        #intrinsic growth rate
                        'R':                        0.5,
                        #intrinsic birth rate (MUST BE 0<=b<=1)
                        'b':                        1,
                        #expectation of distr of n offspring per mating pair
                        #NOTE: using a lambda of equal to
                        # avg_clutch_size * avg_n_clutches_per_yr * avg_surv_rate_to_yr_1
                        # which I put at 4.464 * 2 * 0.16 = 1.428
                        # (where the first value is avg of all CA pops in Tinkle et al 1993,
                        # second value is from Tinkle 1973 and Tinkel et al 1993, and
                        # third comes from the Mt. Diablo population in Ruth 1978,
                        #cited in Tinkle et al. 1993)
                        'n_births_distr_lambda':    1.428,
                        #whether n births should be fixed at n_births_dist_lambda
                        'n_births_fixed':           False,
                        #radius of mate-search area
                        #NOTE: just have individuals mate close by (i.e. within
                        #12.457m, the average interannual movement distance
                        'mating_radius':            0.1
                        }, # <END> 'mating'

                ##########################################
                #### spp num. 0: mortality parameters ####
                ##########################################

                    'mortality'     : {
                        #maximum age
                        # NOTE: 8 is the max age mentioned in Stebbins 1948's
                        # Lassen population study
                        'max_age':                      8,
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
                        #NOTE: the movement and dispersal Wald distributions
                        # are parameterized based on the finding of
                        # Stebbins 1948 that these lizards have extremely
                        # small home ranges;
                        # expressed in cell-widths (where the cell width is
                        # ~730.984 m, according to raster res converted
                        # from dec. deg to m; see NOTE at top of file);
                        # so becuase the mean dispersal distance of
                        # Stebbins' lizards in that study is 12.457m (where
                        # "within/at edge of home area was treated as 0),
                        # this would be 0.01704 cell widths;
                        # then I assigned a sigma that allows for occasional
                        # longer-distance movement events,
                        # but not for long-distance offspring dispersal events
                        #whether or not the species is mobile
                        'move':                     True,
                        #mode of distr of movement direction
                        'direction_distr_mu':       1,
                        #concentration of distr of movement direction
                        'direction_distr_kappa':    0,
                        #mean of distr of movement distance
                        'distance_distr_mu':        0.01704,
                        #variance of distr of movement distance
                        'distance_distr_sigma':     0.01704 / 100,
                        #mean of distr of dispersal distance
                        'dispersal_distr_mu':       0.01704,
                        #variance of distr of dispersal distance
                        'dispersal_distr_sigma':    0.01704 * 100,
                        'move_surf'     : {
                            #move-surf Layer name
                            'layer':                'hab',
                            #whether to use mixture distrs
                            'mixture':              True,
                            #concentration of distrs
                            'vm_distr_kappa':       12,
                            #length of approximation vectors for distrs
                            'approx_len':           5000,
                            }, # <END> 'move_surf'

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
                        #whether starting allele frequencies should be fixed
                        #at 0.5
                        'start_p_fixed':            True,
                        #genome-wide per-base neutral mut rate (0 to disable)
                        'mu_neut':                  0,
                        #genome-wide per-base deleterious mut rate
                        #(0 to disable)
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
                        'allow_ad_hoc_recomb':      False,
                        #whether to save mutation logs
                        'mut_log':                  False,

                        'traits': {

                            ###########################
                            ####trait 0 parameters ####
                            ###########################
                            #trait name (TRAIT NAMES MUST BE UNIQUE!)
                            'thermal_tol': {
                                #trait-selection Layer name
                                'layer':                'tmp',
                                #polygenic selection coefficient
                                'phi':                  0.1,
                                #number of loci underlying trait
                                'n_loci':               10,
                                #mutation rate at loci underlying trait
                                'mu':                   0,
                                #mean of distr of effect sizes
                                'alpha_distr_mu' :      0.1,
                                #variance of distr of effect size
                                'alpha_distr_sigma':    0,
                                #max allowed magnitude for an alpha value
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


    ##########################################################################

    ###############
    #### MODEL ####
    ###############
        'model': {
            #total Model runtime (in timesteps)
            'T':            600,
            #min burn-in runtime (in timesteps)
            'burn_T':       50,
            #seed number
            'num':          None,

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


            ####################################
            #### data-collection parameters ####
            ####################################
            'data': {
                'sampling': {
                    #sampling scheme {'all', 'random', 'point', 'transect'}
                    'scheme':               'transect',
                    #sample size at each point, for point & transect sampling
                    'n':                    250,
                    #coords of collection points, for point sampling
                    'points':               None,
                    #coords of transect endpoints, for transect sampling
                    'transect_endpoints':   [(25, 25), (75, 5)],
                    #num points along transect, for transect sampling
                    'n_transect_points':    6,
                    #collection radius around points, for point & transect
                    #sampling
                    'radius':               10,
                    #when to collect data
                    'when':                 100,
                    #whether to save current Layers when data is collected
                    'include_landscape':    False,
                    #whether to include fixed loci in VCF files
                    'include_fixed_sites':  True,
                    },
                'format': {
                    #format for genetic data {'vcf', 'fasta'}
                    'gen_format':           ['vcf', 'fasta'],
                    #format for vector geodata {'csv', 'shapefile', 'geojson'}
                    'geo_vect_format':      'csv',
                    #format for raster geodata {'geotiff', 'txt'}
                    'geo_rast_format':      'geotiff',
                    },
                }, #<END> 'data'


            } # <END> 'model'

        } # <END> params

    return params


def _run(params, save_figs=False, time_it=False, make_gif=False):
    # set the amount of time before and after climate change
    t_before_cc = 500
    t_after_cc = 100

    # create colormap to match the phenotype colors
    z_cmap = mpl.cm.coolwarm

    # start timer
    if time_it:
        start = time.time()

    # create the model
    from .. import make_model
    mod = make_model(params, verbose=True)

    # set plotting params
    ms = 6
    ax_fontdict = {'fontsize': 12,
                   'name': 'Bitstream Vera Sans'}
    ttl_fontdict = {'fontsize': 15,
                    'name': 'Bitstream Vera Sans'}

    # burn in, then plot starting population, on both rasters
    mod.walk(20000, 'burn')

    # run for the first 500 timsteps, before climate change
    for _ in range(t_before_cc):
        mod.walk(1)
        if save_figs and make_gif:
            save_gif_img(mod.t)

    # calculate neigh-mean raster before climate-change
    neigh_mean_b4 = calc_neighborhood_mean_phenotype(mod)

    # walk for 100 more timesteps, then plot again,
    # at end of climate-change period
    for _ in range(t_after_cc):
        mod.walk(1)
        if save_figs and make_gif:
            save_gif_img(mod.t)

    # end timer
    if time_it:
        stop = time.time()
        tot_time = stop-start

    # calculate neigh-mean raster after climate-change
    neigh_mean_af = calc_neighborhood_mean_phenotype(mod)

    # calculate the difference between the two neigh-mean rasters
    neigh_mean_diff = neigh_mean_af - neigh_mean_b4

    # create and save a population-size plot
    mod.walk(50)
    Nt_fig = plt.figure()
    Nt_ax = Nt_fig.add_subplot(111)
    burn_len = mod.burn_t
    line_height = int(10000*np.ceil(max(mod.comm[0].Nt)/10000))
    Nt_ax.plot(range(len(mod.comm[0].Nt)), mod.comm[0].Nt)
    Nt_ax.plot([burn_len, burn_len], [0, line_height], c='red')
    Nt_ax.plot([burn_len+500, burn_len+500], [0, line_height], c='red')
    Nt_ax.plot([burn_len+600, burn_len+600], [0, line_height], c='red')
    chng_lyr_path = os.path.join(DATA_PATH, 'yosemite_lyrs/ppt/')
    chng_yrs = [int(x[:3]) for x in os.listdir(chng_lyr_path)]
    for yr in chng_yrs:
        Nt_ax.plot([burn_len+yr, burn_len+yr], [0, line_height], ':r',
                   linewidth=0.5)
    Nt_ax.set_label('time (time steps/years)')
    Nt_ax.set_ylabel('total population size (individuals)')
    if save_figs:
        plt.savefig('yosemite_time_series_pop_size.png', format='png', dpi=1000)

    # create the change-raster plot and the draped-raster plot
    make_change_plt(mod)
    DEM = rio.open(os.path.join(DATA_PATH, 'yosemite_DEM.tif')).read()[0, :, :]
    drape_fig = drape_raster(mod, neigh_mean_b4, DEM, cmap='rainbow')

    # create the GIF using imagemagick
    if make_gif:
        try:
            os.system('cd %s' % gif_dir)
            os.system(('cd %s; convert -delay 5 -loop 0 `ls -v` '
                       '../yosemite.gif; cd ..') % gif_dir)
            os.system('..')
        except Exception as e:
            print(('\nCould not use Imagemagick to create the GIF. '
                   'The following error was thrown:\n\t%s') % e)


    # print out time
    if time_it:
        print("\n\nModel ran in %0.2f seconds." % tot_time)

    return mod, drape_fig, neigh_mean_diff, neigh_mean_b4
