#!/usr/bin/python

import numpy as np

params = {

'set_seed' : True,                  #set the seed (for reproducibility)?

'seed_num' : 2,                    #number to seed random number generators

'T' : 2500,                      #total model runtime

'burn_T_min': 50,                     #total burn-in runtime

'L' : 1000,                         #total number of loci

'l_c': [500,250,250],     #chromosome lengths

'recomb_array' : None,           #predetermined linkage map to use for the simulated loci (optional)

'x' : 2,                         #ploidy (for now, leave at 2 for diploidy)

'mu' : 10e-9,                    #genome-wide mutation rate, per base per generation

'alpha_r' : 0.5,                #alpha for beta distribution of linkage values  #NOTE: alpha = 14.999e9, beta = 15e9 has a VERY sharp peak on D = 0.4998333, with no values exceeding equalling or exceeding 0.5 in 10e6 draws in R

'beta_r' : 400,                 #beta for beta distribution of linkage values

'alpha_mut_s' : 25,                # alpha param for the beta distribution describing the highly advantageous selection coeffs for mutations

'beta_mut_s' : 0.5,                # beta param for the beta distribution describing the highly advantageous selection coeffs for mutations

'use_dom' : False,              #whether or not to use dominance (default to False)
                                #NOTE: REALLY JUST NEED TO GET RID OF THE DOMINANCE THING; IT'S ALL MESSED UP

'pleiotropy' : True,           #allow pleiotropy? (i.e. allow same locus to affect value of more than one trait?) false

'traits' : {'num': 2,           #number of traits to simulate
            'scape_num': [0,0], #list of the landscape numbers to be used for selection on each trait 
                                #(i.e.  list length should equal value of 'num' on previous line, as should lengths of subsequent lists)

            'phi' : [0.1, 0.1],   #phenotypic selection coefficient for trait; for each trait, can either be a numeric value, or can be an array of spatially-contingent s values of dimensions equal to land.dims
            'n': [1,100], #number of loci assigned to trait
            'alpha_dist_sigma': [0.5, 0.00025],  #NOTE: for sigma = 0.5, one average locus is enough to generate an both optimum genotypes; for 0.025, 10 loci should (on average, but depends of course on the random sample of alphas drawn!); and so on linearly
            'fitness_fn_gamma': [1, 1], #gamma exponent for the trait's fitness function (<1 = concave up, 1 = linear, >1 = convex up)

            'univ_advant' : [False, False] #is the phenotype unviersally advantageous? if so, phenotypes closer to 1 will have higher fitness at all locations
            },

'N' : 100,                        #total pop size

'dims' : (10,10),             #dimensions of landscape  

'num_scapes' : 2,               #number of landscapes desired

'rand_land' : True,        #whether or not to generate random landscapes
#'rand_land' : False,

'n_rand_pts' : 2600,           #number of random coordinates to be used in generating random landscapes (only needed if rand_land = True)

'islands' : True,           #create habitat islands (outside which individuals cannot move without dying)?

'island_val' : 0.1,              #if greater than 0 (and of course less than 1), the value will be used to
                               #create habitat islands in a random landscape (i.e. all cells less than this 
                               #value will be made completely uninhabitable)

'island_mask' : None,         #pointer to a pickled numpy array, where 0s mask cells outside viable habitat (i.e. outside 'islands'); to use, provide filepath here and set params['island_val'] = 0


'landscape_pt_coords': np.array([[0,0], [0,100], [100,0], [50,40], [100,100], [30,0], [0,30], [70,100], [100,70]]),
#coords of points to use to interpolate defined landscape layers (can be provided as either a single nx2 Numpy array, where n matches the number of points in landscape_pt_vals arrays, to be used as the points for each landscape layer, or a list or tuple of nx2 Numpy arrays, one for each landscape layer; only needed if rand_land = False)

'landscape_pt_vals': [np.array([-0.5,0.5,0.5,0.5,1.5, -0.5, -0.5, 1.5, 1.5]), np.array([1.5,-0.5,-0.5,-0.5,1.5, 0.85, 0.85, 0.85, 0.85])],
#point values to use to interpolate defined landscape layers (a list or tuple of 1xn Numpy arrays, where n matches the number of points in landscape_pt_coords arrays; only needed in rand_land = False)


#'interp_method' : ['nearest'],
'interp_method' : ['linear', 'linear'],   # list of interpolation methods for generation of random landscapes, 1 per landscape to be generated (as set by num_scapes)

'K_cap' : 0.4,                        #per-cell highest carrying capacity value to be reached during burn-in

'move' : True,                     #is this a mobile species?

'movement_surf' : True,       #use a landscape layer as a resistance surface, or habitat quality layer, to direct movement?
#'movement_surf' : False,

'n_movement_surf_scape' : 1,               #scape number to use as the movement surface

'movement_surf_vonmises_kappa' : 2, #kappa value to use in the von Mises mixture distributions (KDEs) underlying resistance surface movement

'movement_surf_gauss_KDE_bandwidth' : 0.2, #bandwidth value to use in the Gaussian KDEs that are created to approximate the von Mises mixture distributions (KDEs) underlying resistance surface movement

'mu_direction' : 0,                #mu for von mises distribution defining movement directions

'kappa_direction' : 0,             #kappa for von mises distribution

'mu_distance' : 0.5,               #mean movement-distance (lognormal distribution)

'sigma_distance' : 0.5,            #sd of movement distance

'sex' : False,                      #is this a sexual species? 

'repro_age' : 0,          #age at sexual maturity (int or float for non-sexual species, tuple or list of two ints/floats for sexual species; set to 'None' to not make this an age-structured species

'dist_weighted_birth' : False,    #should the probability of birth be weighted by the distance between individuals in a pair?

'R':0.5,                            #pop intrinsic growth rate

'b' : 0.2,                         #population intrinsic birth rate (implemented as the probability that an identified potential mating pair successfully mates); 
                                   #NOTE: this may later need to be re-implemented to allow for spatial variation in intrinsic rate (e.g. expression as a raster) and/or for density-dependent births as well as deaths

'd_min' : 0.01,                     #minimum neutral (i.e. non-selection driven) probability of death

'd_max' : 0.90,                     #maximum neutral probability of death

'lambda_offspring': 4,               #expected value of offspring for a successful mating pair (used as the lambda value in a Poisson distribution)

'mating_radius' : 0.5,              #radius of mate-searching area

'mu_dispersal' : 0.5,           #mean dispersal distance (lognormal distribution)

'sigma_dispersal' : 0.5,          #sd of dispersal distance

'size' : 1,              # float/int, or list/tuple of length T containing floats/ints, expressing the target population size over model time as a ratio of the starting size (N)

'sigma_deaths' : 0.2,              # std for the normal distribution used to choose the r.v. of deaths per timestep (mean of this distribution is the overshoot, as calculated from pop.size and pop.census())

'density_dependent_fitness' : True, #should fitness be density dependent? (note: helps to avoid subpopulation 'clumping')


'custom_fns' : {'recomb_rate_fn': None  #if provided, must be a function that returns a single recombination rate value (r) when called
            
               },







'data' : {                      #dictionary defining the data to be collected, the sampling strategy to use, the timesteps for collection, and other parameters

        'sampling_scheme' : 'random',  # can be 'all', 'random', 'point', or 'transect'
        'sampling_args' : {'n':50
                            },
        
                                #args to be unpacked into sampling function (see docstring of sample_data function in data module for details)
        'freq': 15,               #can be an integer (in which case data will be collected every that many timesteps, plus at the end)
                                  #or a list of specific timesteps

        'include_land': False,    #if True, will save the Landscape_Stack object each time other data is saved (probably only useful if land is changing in some way not manually coded by the user)

        'gen_data_format': 'VCF',  #can be 'VCF', 'FASTA', or 'ms'

        'geo_data_format': ['CSV', 'Geotiff'],
                                 #1st argument for points, 2nd for raster; 
                                 #currently 1.) CSV, Shapefile and 2.) Geotiff available

        'run_name': 'test',                #a name for this parameterization and run (used to name the data_directory and files)
        'write_intermittent': True,
        'drop_after_write': True
        },






'stats' : {                      #dictionary defining which stats to be calculated, and parameters on their calculation (including frequency, in timesteps, of collection)
                                 #valid stats include:
                                    # 'Nt'  : population census size
                                    # 'het' : heterozygosity
                                    # 'MAF' : minor allele frequency
                                    # 'LD'  : linkage disequilibrium
        'Nt' : {'calc' : True, 'freq': 1},

        'het' : {'calc' : True, 'freq': 1},

        'MAF' : {'calc' : True, 'freq': 1},
        
        'LD'  : {'calc' : True, 'freq': 1}

        }


}
