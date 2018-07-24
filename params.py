#!/usr/bin/python

import numpy as np


params = {

############
### SEED ###
############

    'seed' : {
        'set_seed'      : True,                     #set the seed? (for reproducibility)
        'seed_num'      : 23432                     #value used to seed random number generators
        }, # <END> 'seed'


###############
#### MODEL ####
###############

    'model' : {
        'burn' : {
            'burn_T_min'    : 50                        #minimum burn-in runtime
            },
        'main' : {
            'T'             : 2500                      #total model runtime
            }
        }, # <END> 'model'


##############
#### LAND ####
##############

    'land' : {
        'dim'          : (25,25),                   #x- and y-dimensionality of landscape  
        'res'          : (1,1),                     #landscape resolution in x and y dimensions (for crosswalking with real-world distances; defaults to meaningless (1,1), will be reset if GIS rasters are read in 
        'ulc'          : (0,0),                     #upper-left corner of the landscape; defaults to meaningless (0,0), can be set to alternative values for crosswalking with real-world data, and will be reset if GIS rasters are read in 
        'n_scapes'    : 2,                        #number of landscapes desired
        'rand_land'     : True,                     #whether or not to generate random landscapes
        'interp_method' : ['linear', 'linear'],     # list of interpolation methods for generation of random landscapes, 
        'n_rand_pts'    : 50,                     #number of random coordinates to be used in generating random landscapes 
                                                        #(only needed if rand_land = True)
                                                        #1 per landscape to be generated (as set by n_scapes)
        'landscape_pt_coords': np.array([[0,0], [0,100], [100,0], [50,40], [100,100], [30,0], [0,30], [70,100], [100,70]]),
                                                    #coords of points to use to interpolate defined landscape layers (can be provided as 
                                                        #either a single nx2 Numpy array, where n matches the number of points in 
                                                        #landscape_pt_vals arrays, to be used as the points for each landscape layer, 
                                                        #or a list or tuple of nx2 Numpy arrays, one for each landscape layer; 
                                                        #only needed if rand_land = False)

        'landscape_pt_vals': [np.array([-0.5,0.5,0.5,0.5,1.5, -0.5, -0.5, 1.5, 1.5]), np.array([1.5,-0.5,-0.5,-0.5,1.5, 0.85, 0.85, 0.85, 0.85])],
                                                    #point values to use to interpolate defined landscape layers (a list or tuple of 
                                                        #1xn Numpy arrays, where n matches the number of points in landscape_pt_coords arrays;
                                                        #only needed if rand_land = False)
        'density_grid_window_width': None,          #with window-width used for the Density_Grid_Stack used to calculate
                                                        #population density (if set to None, defaults to the closest factor 
                                                        #of the larger landscape dimension to 1/10th of that dimension)
        'gis'           : {
                    'scape_nums'                 : [1],   #list of the scape_nums for which GIS layers should be read in
                    'filepaths'                  : [], #'/home/ihavehands/Desktop/stuff/berk/research/projects/sim/yos_30yr_normals_90x90.tif'], #list of the filepaths to read into the scapes indicated by the scape_num list
                    'scale_min_vals'                   : [-1.37],    #minimum values to use for rescaling the rasters (will be rescaled to 0<=x<=1); may be different than the actual minimum values in the rasters, if they will be changing to future rasters with values outside the range of these rasters
                    'scale_max_vals'                   : [19.11]   #maxmimum input values against which to rescale the rasters to (will be rescaled to 0<=x<=1)
                    },
        'islands'       : {
                    'islands'                   : False,                    #create habitat islands (outside which individuals cannot move without dying)?
                    'island_val'                : 0,                            #if greater than 0 (and of course less than 1), the value will be used to
                                                                                #create habitat islands in a random landscape (i.e. all cells less than this 
                                                                                #value will be made completely uninhabitable)
                    'island_mask'               : None,                     #pointer to a pickled numpy array, where 0s mask cells outside viable habitat 
                                                                                #(i.e. outside 'islands'); to use, provide filepath here and set 
                                                                                #params['island_val'] = 0
                        }, # <END> 'islands'
        'move_surf'     : {
                'move_surf'                 : True,                     #use a landscape layer as a resistance surface, or habitat quality 
                                                                                #layer, to direct movement?
                'move_surf_scape_num'         : 0,                        #scape number to use as the movement surface
                'move_surf_approx_len'      : 7500,                     #length of the lookup vectors (numpy arrays) used to approximate 
                                                                                #the VonMises mixture distributions at each cell
                'move_surf_vonmises_kappa'  : 2,                        #kappa value to use in the von Mises mixture distributions (KDEs) 
                                                                                #underlying resistance surface movement
                'move_surf_gauss_KDE_bandwidth' : 0.2                   #bandwidth value to use in the Gaussian KDEs that are created to 
                                                                                #approximate the von Mises mixture distributions (KDEs) 
                                                                                #underlying resistance surface movement
                            } # <END> 'move_surf'
        }, # <END> 'landscape'
 

################
#### GENOME ####
################

    'genome' : {
        'L'             : 1000,                     #total number of loci
        'l_c'           : [500, 500],               #chromosome lengths [sum(l_c) == L enforced]
        'recomb_array'  : None,                     #predetermined linkage map to use for the simulated loci (optional)
        'x'             : 2,                        #ploidy (for now, leave at 2 for diploidy)
        'mu'            : 10e-9,                    #genome-wide mutation rate, per base per generation
        'alpha_r'       : 0.5,                      #alpha for beta distribution of linkage values  
                                                        #NOTE: alpha = 14.999e9, beta = 15e9 has a VERY sharp peak on D = 0.4998333, 
                                                        #with no values exceeding equalling or exceeding 0.5 in 10e6 draws in R
        'beta_r'        : 400,                      #beta for beta distribution of linkage values
        'alpha_mut_s'   : 25,                       # alpha param for the beta distribution describing the highly advantageous selection coeffs for mutations
        'beta_mut_s'    : 0.5,                      # beta param for the beta distribution describing the highly advantageous selection coeffs for mutations

        'use_dom'       : False,                    #whether or not to use dominance (default to False)
                                                        #NOTE: REALLY JUST NEED TO GET RID OF THE DOMINANCE THING; IT'S ALL MESSED UP
        'pleiotropy'    : True,                     #allow pleiotropy? (i.e. allow same locus to affect value of more than one trait?) false
        'recomb_rate_custom_fn': None,              #if provided, must be a function that returns a single recombination rate value (r) when called
        'recomb_lookup_array_size': int(1e3),       #the size of the recombination-path lookup array to have
        #read in at one time (needs to be comfortably larger than the anticipated totaly number of
        #recombination paths to be drawn at once, i.e. than 2 times the anticipated most number of births at once)
        'n_recomb_paths': int(1e4),               #the total number of distinct recombination paths to
        #generate at the outset, to approximate truly free recombination at the recombination rates specified
        #by the genomic architecture (hence the larger the value the less the likelihood of mis-approximation artifacts)
        'traits'        : {
                    'num'       : 1,                        #number of traits to simulate
                    'scape_num' : [1],                  #list of the landscape numbers to be used for selection on each trait 
                                                                #(i.e.  list length should equal value of 'num' on previous line, 
                                                                #as should lengths of subsequent lists)
                    'phi' : [0.1],             #phenotypic selection coefficient for trait; 
                                                                #for each trait, can either be a numeric value, or can be an array 
                                                                #of spatially-contingent s values of dimensions equal to land.dims
                    'n_loci': [10],                  #number of loci assigned to trait
                    'alpha_dist_sigma': [0.5],    #NOTE: for sigma = 0.5, one average locus is enough to generate both optimum 
                                                                #genotypes; for 0.025, 10 loci should (on average, but depends of course on 
                                                                #the random sample of alphas drawn!); and so on linearly
                    'gamma': [1],                     #gamma exponent for the trait's fitness function (<1 = concave up, 1 = linear, >1 = convex up)
                    'univ_advant' : [False]    #is the phenotype unviersally advantageous? 
                                                                #if so, phenotypes closer to 1 will have higher fitness at all locations
                        } # <END> 'traits'
        
    }, # <END> 'genome'


#############
#### POP ####
#############

    'pops' : {
        0  :   {
            'name'       : 'pop0',         #each population can take a string as a name
            'start'      : {
               'N'       : 200,             #starting population size
               'K_scape_num'   : 0,                #the scape_num of the raster to use as the carrying-capacity raster (K)
               'K_fact'        : 2                 #the factor to multiply the K raster by in order to generate pop.K
               }, # <END> 'main'

            'mating'    : {
               'repro_age'     : 0,                    #age at sexual maturity (int or float for non-sexual species, tuple or list 
                                                           #of two ints/floats for sexual species; set to 'None' to not make this 
                                                           #an age-structured species
               'max_age'       : 3,                      #age beyond which all individuals will automatically die; default to None
               'sex'           : False,                #is this a sexual species? 
               'dist_weighted_birth' : False,          #should the probability of birth be weighted by the distance between 
                                                           #individuals in a pair?
               'R'             : 0.5,                  #pop intrinsic growth rate
               'b'             : 0.2,                  #population intrinsic birth rate (implemented as the probability 
                                                           #that an identified potential mating pair successfully mates); 
                                                           #NOTE: this may later need to be re-implemented to allow for spatial 
                                                           #variation in intrinsic rate (e.g. expression as a raster) and/or for 
                                                           #density-dependent births as well as deaths
               'n_births_lambda': 4,                   #expected value of offspring for a successful mating pair (used as the lambda value in a Poisson distribution)
               'mating_radius' : 1                   #radius of mate-searching area
               }, # <END> 'mating'

            'mortality'     : {
               'n_deaths_sigma'            : 0.2,      # std for the normal distribution used to choose the r.v. of deaths 
                                                           #per timestep (mean of this distribution is the overshoot, 
                                                           #as calculated from pop.size and pop.census())
               'density_dependent_fitness' : True,     #should fitness be density dependent? (note: helps to avoid 
               'd_min'                     : 0.01,     #minimum neutral (i.e. non-selection driven) probability of death
               'd_max'                     : 0.90,     #maximum neutral probability of death
               'islands'                   : False

        }   , # <END> 'mortality'
                                                           #subpopulation 'clumping')
            'movement'   : {
               'movement'          : True,                 #is this a mobile species?
               'direction_mu'      : 0,                    #mu for von mises distribution defining movement directions
               'direction_kappa'   : 0,                    #kappa for von mises distribution
               'distance_mu'       : 0.5,                  #mean movement-distance (lognormal distribution)
               'distance_sigma'    : 0.5,                  #sd of movement distance
               'dispersal_mu'      : 0.5,                  #mean dispersal distance (lognormal distribution)
               'dispersal_sigma'   : 0.5,                  #sd of dispersal distance
               'move_surf'     : True                  #use the movement surface for this population?
            }    # <END> 'movement'
        } # <END> '0'

            
    #NOTE: COPY AND PASTE 0 HERE AND GIVE A DIFFERENT NAME, TO CREATE ADDITIONAL POPULATIONS IF DESIRED
   

    }, # <END> 'pops'


################
#### CHANGE ####
################

    'change' : {
            'pops' : {
                    'dem'   :   {                       
                        'pop0'  : {                     #to which population should the following changes apply?
                                                            #(all population sizes are expressed relative to the carrying-capacity 
                                                            #raster at the time that the demographic change event begins (i.e. as
                                                            #factors by which pop.K will be multiplied; thus they can also be thought 
                                                            #of as multipliers of the expected total population size (i.e. pop.K.sum())
                                                            #and they will generally change the average population size by that multiple, 
                                                            #but of course not precisely, because population size is stochastic. If you 
                                                            #seek exact control of total population size, please seek a simpler simulation 
                                                            #model, perhaps a coalescent one.

                            0 : {                               #can add an arbitrary number of demographic change events for 
                                                                    #each population, each event identified by a distinct integer
                                'kind'          : 'custom',  #what kind of change event? ('monotonic', 'stochastic', 'cyclical', 'custom')
                                'start'         : 200,          #at which timestep should the event start?
                                'end'           : 1200,          #at which timestep should the event end?
                                'rate'          : .98,         #at what rate should the population change each timestep 
                                                                    #(positive for growth, negative for reduction)
                                'interval'      : 11,         #at what interval should stochastic change take place (None defaults to every timestep)
                                'dist'          : 'uniform',    #what distribution to draw stochastic population sizes from (valid values: 'uniform', 'normal')
                                'size_target'   : None,         #what is the target size of the demographic change event (defaults to None)
                                'n_cycles'      : 20,         #how many cycles of cyclical change should occur during the event?
                                'size_range'    : (0.5, 1.5),     #a tuple of the min and max population sizes to be used in stochastic or cyclical changes
                                'timesteps'     : [10,100,1000,1100,1150,1200],         #at which timesteps should custom changes take place?
                                'sizes'         : [3,2,1.5,2,0.5,1]          #what custom size-changes should occur at the above-stipulated timesteps?
                                }
                            }
                        },
                    'other' :   {       
                        'pop0'  : {                     #to which population should the following changes apply?
                            'b' : {                         #keys are the parameters to be changed, and values
                                                                #are dictionaries containing a list of timesteps at which to changed those
                                                                #parameters and a list of values to which to change them
                                'timesteps'     : None,
                                'values'        : None
                                }
                            }
                        }
                    },
            'land' : {
                      1 : {
                            'end_rast' : np.zeros((90,90)),
                            't_start' : 1500,
                            't_end'   : 2000,
                            'n'       : 10
                            },
                      0 : {
                            'end_rast' : np.zeros((90,90)),
                            't_start' : 1500,
                            't_end'   : 2000,
                            'n'       : 2
                            }
                     
                    },
            'gen' : {
                }
        }, # <END> 'change'


##############
#### DATA ####
##############

    'data' : {      #dictionary defining the data to be collected, the sampling strategy to use, the timesteps for collection, and other parameters

        'sampling_scheme'   : 'random',             #can be 'all', 'random', 'point', or 'transect'
        'sampling_args'     : {                     #args to be unpacked into sampling function (see docstring of sample_data function in data module for details)
                        'n'     : 50        
                            },
        'freq'              : 15,                   #can be an integer (in which case data will be collected every that many timesteps, plus at the end)
                                                        #or a list of specific timesteps
        'include_land'      : False,                #if True, will save the Land object each time other data is saved (probably only useful 
                                                        #if land is changing in some way not manually coded by the user)
        'gen_data_format'   : 'VCF',                #can be 'VCF', 'FASTA', or 'ms'
        'geo_data_format'   : ['CSV', 'Geotiff'],   #1st argument for points, 2nd for raster; currently 1.) CSV, Shapefile and 2.) Geotiff available
        'run_name'          : 'test',               #a name for this parameterization and run (used to name the data_directory and files)
        'write_intermittent': True,
        'drop_after_write'  : True
        }, #<END> 'data'


###############
#### STATS ####
###############

    'stats' : {     #dictionary defining which stats to be calculated, and parameters on their calculation (including frequency, in timesteps, of collection)
                        #valid stats include:
                            # 'Nt'  : population census size
                            # 'het' : heterozygosity
                            # 'MAF' : minor allele frequency
                            # 'LD'  : linkage disequilibrium
        'Nt'                : {'calc' : True, 'freq': 1},
        'het'               : {'calc' : True, 'freq': 1},
        'MAF'               : {'calc' : True, 'freq': 1},
        'LD'                : {'calc' : True, 'freq': 1}
    }, # <END> 'stats'


###############
#### OTHER ####
###############

    'other' : {
    } #<END> 'other'

} # <END> 'params'

# <END> params.py
