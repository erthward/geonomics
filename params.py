#!/usr/bin/python
#params.py

'''Geonomics parameters file.'''


######################################
#TODO:
 # create a params.py MODULE instaed, with a function for generating a template params dictionary, taking
 # arguments for the number of scapes, number of pops, whether or not to include a genome, whether or not
 # to include land and pop changes, whether or not to include a model-manager section with data and/or stats

 #if I decide to use the recursive function below to create a dynamic-attribute dict class of arbitrary depth,
 #then IMPORTANT to check that none of the following serve as keys:
    # ['clear', 'copy', 'fromkeys', 'get', 'items', 'keys', 'pop', 'popitem', 'setdefault', 'update', 'values']

#if I decide to use the Params_Dict class, then keys CANNOT be numbers! (so come up with a different scheme
#for pops and scapes ... perhaps the key should just be the name that the user wants, and names will be
#rejected iff they clobber one of the dict methods listed above?

#also, if I use the Params_Dict class then I could change all the params keys to use some standardization that
#doesn't appear elsewhere in the package, e.g. ALL-CAPS-AND-HYPHENS?

######################################



import numpy as np

params = {

##############
#### LAND ####
##############

    'land': {

    ##############
    #### main ####
    ##############

        'main': {
            'dim':                      (20,20),
                #x- and y-dimensionality of landscape  
            'res':                      (1,1),
                #landscape resolution in x and y dimensions (for crosswalking with real-world 
                #distances; defaults to meaningless (1,1), will be reset if GIS rasters are read in 
            'ulc':                      (0,0),
                #upper-left corner of the landscape; defaults to meaningless (0,0), can be set 
                #to alternative values for crosswalking with real-world data, and will be 
                #reset if GIS rasters are read in 
            'prj':                      None,
                #projection of the landscape; only applicable if layers are to
                #be read in from a raster file; defaults to None
            }, # <END> 'main'

    ################
    #### scapes ####
    ################

        'scapes': {

            0: {
                #this integer key should be incremented for each successive scape

        ##############
        #### init ####
        ##############

                'init': {
                    #initiating parameters for this scape
                    'name':             'scape0',
                        #each scape can take a unique string name (e.g. 'scape0', '1994', 'mean_T')
                    'rand': {
                        #parameters for making a scape using interpolation from randomly located random values
                        'n_pts':                        500,
                            #number of random coordinates to be used in generating random landscapes 
                                #(only needed if rand == True)
                        'interp_method':                'cubic'
                            # interpolation method (valid: 'linear', 'cubic', 'nearest')
                        }
#                    'defined': {
#                        #parameters for making a scape using interpolation from a defined set of valued points
#                        'pts':                    None,
#                            #coords of points to use to interpolate defined scape (provided as 
#                                #a nx2 Numpy array, where n matches the number of points in 
#                                #the scape_pt_vals array, to be used as the points
#                                #to be interpolated; only needed if rand == False)
#                        'vals':                      None,
#                            #point values to use to interpolate defined landscape layers (a  1xn Numpy array, 
#                                #where n matches the number of points in scape_pt_coords arrays;
#                                #only needed if rand == False)
#                        'interp_method':                None
#                            # interpolation method (valid: 'linear', 'cubic', 'nearest')
#                        },
#                    'gis': {
#                        #parameters for making a scape using a GIS raster file
#                        'filepath':                     '/home/ihavehands/Desktop/stuff/berk/research/projects/sim/yos_30yr_normals_90x90.tif', 
#                            #filepath to read into this scape
#                        'scale_min_val':                -1.37,
#                            #minimum values to use for rescaling the raster (will be
#                                #rescaled to 0<=x<=1); NOTE: this may be different than the actual minimum 
#                                #value in the raster, especially if raster will be changing to a future raster 
#                                #with values outside the range of this one
#                        'scale_max_val':                19.11   
#                            #maxmimum input value against to which to rescale the raster (will be rescaled to 0<=x<=1)
#                        }
                    } # <END> 'init'

        ################
        #### change ####
        ################

#                'change': {
#                    #land-change events for this scape
#                    'end_rast':         np.zeros((90,90)),
#                        #scape to be set as the endpoint of the land-change event
#                    'start_t':          1500,
#                        #timestep on which to start the land-change event
#                    'end_t':            2000,
#                        #timestep on which to end the land-change event
#                    'n_steps':          10
#                        #number of stepwise changes to make between t_start and t_end
#                        }
                }, # <END> scape 0


    #*****
    #NOTE: COPY AND PASTE THE scape SECTION HERE AND GIVE A DIFFERENT INDEX TO CREATE ADDITIONAL SCAPES FOR
    #THIS LAND
    #*****

            1: {
                #this integer key should be incremented for each successive scape

        ##############
        #### init ####
        ##############

                'init': {
                    #initiating parameters for this scape
                    'name':             'scape1',
                        #each scape can take a unique string name (e.g. 'scape0', '1994', 'mean_T')
                    'rand': {
                        #parameters for making a scape using interpolation from randomly located random values
                        'n_pts':                        500,
                            #number of random coordinates to be used in generating random landscapes 
                                #(only needed if rand == True)
                        'interp_method':                'nearest'
                            # interpolation method (valid: 'linear', 'cubic', 'nearest')
                        }
#                    'defined': {
#                        #parameters for making a scape using interpolation from a defined set of valued points
#                        'pts':                    None,
#                            #coords of points to use to interpolate defined scape (provided as 
#                                #a nx2 Numpy array, where n matches the number of points in 
#                                #the scape_pt_vals array, to be used as the points
#                                #to be interpolated; only needed if rand == False)
#                        'vals':                      None,
#                            #point values to use to interpolate defined landscape layers (a  1xn Numpy array, 
#                                #where n matches the number of points in scape_pt_coords arrays;
#                                #only needed if rand == False)
#                        'interp_method':                None
#                            # interpolation method (valid: 'linear', 'cubic', 'nearest')
#                        },
#                    'gis': {
#                        #parameters for making a scape using a GIS raster file
#                        'filepath':                     '/home/ihavehands/Desktop/stuff/berk/research/projects/sim/yos_30yr_normals_90x90.tif', 
#                            #filepath to read into this scape
#                        'scale_min_val':                -1.37,
#                            #minimum values to use for rescaling the raster (will be
#                                #rescaled to 0<=x<=1); NOTE: this may be different than the actual minimum 
#                                #value in the raster, especially if raster will be changing to a future raster 
#                                #with values outside the range of this one
#                        'scale_max_val':                19.11   
#                            #maxmimum input value against to which to rescale the raster (will be rescaled to 0<=x<=1)
#                        }
                    }, # <END> 'init'

        ################
        #### change ####
        ################

                'change': {
                    #land-change events for this scape
                    'end_rast':         np.zeros((20,20)),
                        #scape to be set as the endpoint of the land-change event
                    'start_t':          3,
                        #timestep on which to start the land-change event
                    'end_t':            7,
                        #timestep on which to end the land-change event
                    'n_steps':          2
                        #number of stepwise changes to make between t_start and t_end
                        }
                }, # <END> scape 1

            } # <END> 'scapes' 

        }, # <END> 'land'


###################
#### COMMUNITY ####
###################

    'comm': {

        'pops': {

            0  :   {
                #this integer key should be incremented for each successive population

            ##############
            #### init ####
            ##############

                'init': {
                    'name': 'pop0',
                        #each pop can take a unique string name (e.g. 'pop0', 'south', 'C_fasciata')
                    'N':                200,
                        #starting population size
                    'K_scape_num':      0,
                        #the scape_num of the raster to use as the carrying-capacity raster (K)
                    'K_fact':           2
                        #the factor to multiply the K raster by in order to generate pop.K
                    }, # <END> 'init'

            ################
            #### mating ####
            ################

                'mating'    : {
                    'repro_age':            0,
                        #age at sexual maturity (int or float for non-sexual species, tuple or list 
                            #of two ints/floats for sexual species; set to 'None' to not make this 
                            #an age-structured species
                    'max_age':              5,
                        #age beyond which all individuals will automatically die; default to None
                    'sex':                  False,
                        #is this a sexual species? 
                    'dist_weighted_birth':  False,
                        #should the probability of birth be weighted by the distance between 
                            #individuals in a pair?
                    'R':                    0.5,
                        #pop intrinsic growth rate
                    'b':                    0.2,
                        #population intrinsic birth rate (implemented as the probability 
                            #that an identified potential mating pair successfully mates); 
                            #NOTE: this may later need to be re-implemented to allow for spatial 
                            #variation in intrinsic rate (e.g. expression as a raster) and/or for 
                            #density-dependent births as well as deaths
                    'n_births_lambda':      4,
                        #expected value of offspring for a successful mating pair (used as the lambda value in a Poisson distribution)
                    'mating_radius':        1
                        #radius of mate-searching area
                    }, # <END> 'mating'

            ###################
            #### mortality ####
            ###################

                'mortality'     : {
                    'n_deaths_sigma':           0.2,
                        #std for the normal distribution used to choose the r.v. of deaths 
                            #per timestep (mean of this distribution is the overshoot, 
                            #as calculated from pop.size and pop.census())
                    'selection':                True,
                        #should the population undergo natural selection?
                    'dens_dependent_fitness':   True,
                        #should fitness be density dependent? (note: helps to avoid subpopulation 'clumping')
                    'dens_grid_window_width':   None,
                        #with window-width used for the Density_Grid_Stack that calculates pop density 
                            #(if set to None, defaults to the closest factor of the larger landscape 
                            #dimension to 1/10th of that dimension)
                            #NOTE: will eventually default to an approximation of Wright's genetic neighborhood 
                            #distance, based on the population's movement/dispersal parameters
                   'd_min':                     0.01,
                        #minimum neutral (i.e. non-selection driven) probability of death
                    'd_max':                    0.90,
                        #maximum neutral probability of death
                    'islands':  {
                        'make':                 False,
                            #create habitat islands (outside which individuals cannot move without dying)?
                        'island_val':           0
                            #if greater than 0 (and of course less than 1), the value will be used to
                                #create habitat islands in a random landscape (i.e. all cells less than this 
                                #value will be made completely uninhabitable)
                        }, # <END> 'islands'
 
                    }, # <END> 'mortality'

            ##################
            #### movement ####
            ##################

                'movement': {
                   'move':          True,                 
                        #is this a mobile species?
                    'direction_mu':     0,
                        #mu for von mises distribution defining movement directions
                    'direction_kappa':  0,
                        #kappa for von mises distribution
                    'distance_mu':      0.5,
                        #mean movement-distance (lognormal distribution)
                    'distance_sigma':   0.5,
                        #sd of movement distance
                    'dispersal_mu':     0.5,
                        #mean dispersal distance (lognormal distribution)
                    'dispersal_sigma':  0.5,
                        #sd of dispersal distance
                    'move_surf'     : {
                        'make':                         True,                     
                            #use a landscape layer as a resistance surface, or habitat quality layer, to direct movement?
                        'scape_num':                    0,
                            #scape number to use as the movement surface
                        'approximation_len':            7500,
                            #length of the lookup vectors (numpy arrays) used to approximate
                                #the VonMises mixture distributions at each cell
                        'vm_kappa':                     None,
                            #kappa value to use in the von Mises mixture distributions (KDEs) 
                                #underlying resistance surface movement
                        'gauss_KDE_bw':                 None
                            #bandwidth value to use in the Gaussian KDEs that are created to 
                                #approximate the von Mises mixture distributions (KDEs) 
                                #underlying resistance surface movement
                        } # <END> 'move_surf'

                    },    # <END> 'movement'

            ################
            #### genome ####
            ################

                'genome': {
                    'L':                        1000,
                        #total number of loci
                    'l_c':                      [500,300,200],
                        #chromosome lengths [sum(l_c) == L enforced]
                    'gen_arch_file':            None,
                        #if not None, should point to a file stipulating a
                            #custom genomic architecture (i.e. a CSV with loci 
                            #as rows and 'locus_num', 'p', 'r', 'trait', and 
                            #'alpha' as columns, such as is created by
                            #main.make_params_file, when the custom_gen_arch
                            #arugment is True)
                    'mu_neut':          1e-9,
                        #genome-wide neutral mutation rate, per base per generation
                            #(set to 0 to disable neutral mutation)
                    'mu_delet':            1e-9,
                        #genome-wide deleterious mutation rate, per base per generation
                            #(set to 0 to disable deleterious mutation)
                            #NOTE: these mutations will fall outside the loci involved in any traits
                            #being simulated, and are simply treated as universally deleterious, with the same
                            #negative influence on fitness regardless of spatial context
                    'mut_log':              False,
                        #whether or not to store a mutation log; if true, will be saved as mut_log.txt
                        #within each iteration's subdirectory
                    'shape_delet_s_dist':      0.2,
                    'scale_delet_s_dist':   0.2,
                        #mean and standard deviation of the per-allele effect size of deleterious mutations (std = 0 will fix all
                            #mutations for the mean value)
                    'alpha_r_dist':                  0.5,
                        #alpha for beta distribution of linkage values  
                            #NOTE: alpha = 14.999e9, beta = 15e9 has a VERY sharp peak on D = 0.4998333, 
                            #with no values exceeding equalling or exceeding 0.5 in 10e6 draws in R
                    'beta_r_dist':                   400,
                        #beta for beta distribution of linkage values
                    'use_dom':                  False,
                        #whether or not to use dominance (default to False)
                            #NOTE: REALLY JUST NEED TO GET RID OF THE DOMINANCE THING; IT'S ALL MESSED UP
                    'pleiotropy':               True,
                        #allow pleiotropy? (i.e. allow same locus to affect value of more than one trait?) false
                    'recomb_rate_custom_fn':    None,
                        #if provided, must be a function that returns a single recombination rate value (r) when called
                    'recomb_lookup_array_size': int(1e3),
                        #the size of the recombination-path lookup array to have
                            #read in at one time (needs to be comfortably larger than the anticipated totaly number of
                            #recombination paths to be drawn at once, i.e. than 2 times the anticipated most number of births at once)
                    'n_recomb_paths':           int(1e4),
                        #the total number of distinct recombination paths to
                            #generate at the outset, to approximate truly free recombination at the recombination rates specified
                            #by the genomic architecture (hence the larger the value the less the likelihood of mis-approximation artifacts)

                    'traits': {
                        0: {
                            #an arbitrary number of traits can be provided for a genomic_architecture object
                            'name':             'trait0',
                                #each trait must be a given a string name (e.g. 'trait0', 'scape0_trait', 'min_temp_trait', 'bill_length')
                            'scape_num':        1,
                                #the landscape numbers to be used for selection on this trait 
                            'phi':              0.1,
                                #phenotypic selection coefficient for this trait; can either be a 
                                    #numeric value, or can be an array of spatialized selection 
                                    #values (with dimensions equal to land.dims)
                            'n_loci':           1,
                                #number of loci to be assigned to this trait
                            'mu':      1e-9,
                                #mutation rate for this trait (if set to 0, or if genome['mutation'] == False, no mutation will occur)
                                    #(set to 0 to disable mutation for this trait)
                            'mean_alpha_dist' : 0,
                            'std_alpha_dist' : 0.5,
                                #the mean and standard deviation of the normal distribution used to choose effect size
                                    #(alpha) for this trait's loci
                                    #NOTE: for mean = 0, std = 0.5, one average locus is enough to generate both optimum 
                                    #genotypes; for mean = 0, std = 0.025, 10 loci should generate both (on average, but depends of course on 
                                    #the random sample of alphas drawn); and so on linearly
                            'gamma':            1,
                                #gamma exponent for the trait's fitness function (determines the shape of the
                                #curve of fitness as a function of absolute difference between an individual's
                                #phenotype and its environment; <1 = concave up, 1 = linear, >1 = convex up)
                            'univ_advant':      False
                                #is the trait unviersally advantageous? if so, phenotypes closer to 1 will 
                                    #have higher fitness at all locations on the land
                            }, # <END> trait 0

                        
                        1: {
                            #an arbitrary number of traits can be provided for a genomic_architecture object
                            'name':             'trait1',
                                #each trait must be a given a string name (e.g. 'trait0', 'scape0_trait', 'min_temp_trait', 'bill_length')
                            'scape_num':        1,
                                #the landscape numbers to be used for selection on this trait 
                            'phi':              0.1,
                                #phenotypic selection coefficient for this trait; can either be a 
                                    #numeric value, or can be an array of spatialized selection 
                                    #values (with dimensions equal to land.dims)
                            'n_loci':           10,
                                #number of loci to be assigned to this trait
                            'mu':      1e-9,
                                #mutation rate for this trait (if set to 0, or if genome['mutation'] == False, no mutation will occur)
                            'mean_alpha_dist' : 0,
                            'std_alpha_dist' : 0.5,
                                #the mean and standard deviation of the normal distribution used to choose effect size
                                    #(alpha) for this trait's loci
                                    #NOTE: for mean = 0, std = 0.5, one average locus is enough to generate both optimum 
                                    #genotypes; for mean = 0, std = 0.025, 10 loci should generate both (on average, but depends of course on 
                                    #the random sample of alphas drawn); and so on linearly
                            'gamma':            1,
                                #gamma exponent for the trait's fitness function (determines the shape of the
                                #curve of fitness as a function of absolute difference between an individual's
                                #phenotype and its environment; <1 = concave up, 1 = linear, >1 = convex up)
                            'univ_advant':      False
                                #is the trait unviersally advantageous? if so, phenotypes closer to 1 will 
                                    #have higher fitness at all locations on the land
                            }, # <END> trait 1

        #*****
        #NOTE: COPY AND PASTE THE trait SECTION HERE AND GIVE A DIFFERENT NAME TO CREATE ADDITIONAL TRAITS
        #FOR THIS GENOMIC ARCHITECTURE
        #*****

                        }, # <END> 'traits'

                    }, # <END> 'genome'

            ################
            #### change ####
            ################

                'change': {
                    'dem': {                       
                        #(all population sizes are expressed relative to the carrying-capacity 
                            #raster at the time that the demographic change event begins (i.e. as
                            #factors by which pop.K will be multiplied; thus they can also be thought 
                            #of as multipliers of the expected total population size (i.e. pop.K.sum())
                            #and they will generally change the average population size by that multiple, 
                            #but of course not precisely, because population size is stochastic. If you 
                            #seek exact control of total population size, please seek a simpler simulation 
                            #model, perhaps a coalescent one.
                        0: {
                            #can add an arbitrary number of demographic change events for 
                                #each population, each event identified by a distinct integer
                            'kind':             'custom',
                                #what kind of change event? ('monotonic', 'stochastic', 'cyclical', 'custom')
                            'start':            200,
                                #at which timestep should the event start?
                            'end':              1200,
                                #at which timestep should the event end?
                            'rate':             .98,
                                #at what rate should the population change each timestep 
                                    #(positive for growth, negative for reduction)
                            'interval':         11,
                                #at what interval should stochastic change take place (None defaults to every timestep)
                            'dist':             'uniform',
                                #what distribution to draw stochastic population sizes from (valid values: 'uniform', 'normal')
                            'size_target':      None,
                                #what is the target size of the demographic change event (defaults to None)
                            'n_cycles':         20,
                                #how many cycles of cyclical change should occur during the event?
                            'size_range':       (0.5, 1.5),     
                                #an iterable of the min and max population sizes to be used in stochastic or cyclical changes
                            'timesteps':        [6,8],
                                #at which timesteps should custom changes take place?
                            'sizes':            [2,0.25],
                                #what custom size-changes should occur at the above-stipulated timesteps?
                            } # <END> event 0
                            
                        }, # <END> 'dem'

                    'parameters': {       
                        #other (i.e. non-demographic) population change events
                        'b': {
                            #the parameters to be changed should be the keys in this dict, 
                                #and values are dictionaries containing a list of timesteps 
                                #at which to changed those parameters and a list of values 
                                #to which to change them
                            'timesteps':        None,
                            'vals':           None
                                }
                            } # <END> 'other'

                        } # <END> 'change'

                    }, # <END> pop 0

        #*****
        #NOTE: COPY AND PASTE THE pop SECTION HERE AND GIVE A DIFFERENT NAME TO CREATE ADDITIONAL
        #POPULATIONS FOR THIS COMMUNITY
        #*****

            }, # <END> 'pops'

        }, # <END> 'comm'


###############
#### MODEL ####
###############

    'model': {
        'seed': {
            #parameters to control whether and how to set the seed
            'set':          True,
                #set the seed? (for reproducibility)
            'num':          304
                #value used to seed random number generators
            }, # <END> 'seed'
        'its': {
            #parameters to control how many iterations of the model to run,
            #and whether or not to randomize the land and/or community 
            #objects in each model iteration
            'n_its': 5,
                #how many iterations of the model should be run?
            'rand_land':    False,
                #randomize the land for each new iteration?
            'rand_comm':    False,
                #randomize the community for each new iteration?
            'rand_burn':  False,
                #randomize the burn-in for each new iteration? (i.e. burn in
                #each time, or burn in once at creation and then use the same
                #burnt-in population for each iteration?)
            }, # <END> 'iterations'
        'time': {
            #parameters to control the number of burn-in and main timesteps to
            #run for each iterations
            'T':            40,
                #total model runtime (in timesteps)
            'burn_T':       30
                #minimum burn-in runtime (in timesteps; this is a mininimum because 
                    #burn-in will run for at least this long but until
                    #stationarity detected, which will likely be longer)
            }, # <END> 'timesteps'
        'data': {
            #dictionary defining the data to be collected, the sampling 
            #strategy to use, the timesteps for collection, and other parameters
            'sampling': {
                #args to be unpacked into sampling function (see docstring 
                    #of sample_data function in data module for details)
                'scheme':               'random',
                    #valid: 'all', 'random', 'point', or 'transect'
                'n':                    50,
                    #size of samples to be collected (in number of individuals)
                'points':               None,
                    #the x,y points at which data should be sampled (expressed
                        #as a list or tuple of length-2 lists or 2-tuples)
                'transect_endpoints':   [(0,0), (20,20)],
                    #endpoints of the transect to be sampled (only needed if
                        #scheme is 'transect), expressed as a pair of 
                        #ordered x,y pairs (in tuples or lists)
                'n_transect_points':    4,
                    #the number of evenly spaced points along the transect
                    #at which to sample (only needed if scheme is 'transect')
                'radius':               2,
                    #radius around sampling points within which to sample
                    #individuals (only needed is scheme is 'point' or
                    #'transect')
                'when':                 15,
                    #can be an integer (in which case data will be collected every 
                    #that many timesteps, plus at the end) or a list of specific
                    #timesteps; a value of 0 or None will default to a single
                    #data-collection step after the model has run
                'include_land':         True,
                    #if True, will save the Land object each time other data is saved 
                    #(probably only useful if land is changing in some way not manually coded by the user)
                'include_fixed_sites':  False,
                    #if True, and if genetic data is to be formatted as VCFs,
                        #the VCFs will contain fixed sites, not just variants
                        #(defaults to False)
                },
            'format': {
                'gen_format':           ['vcf', 'fasta'],
                    #format to use for saving genetic data; 
                        #currently valid values: 'vcf', 'fasta', 
                        #or a list containing both, if both 
                        #should be written
                'geo_vect_format':      'csv',
                    #format to use for saving geographic points; 
                        #currently valid values: 'csv', 'shapefile', 'geojson'
                'geo_rast_format':      'geotiff',
                    #format to use for saving landscape rasters (which will
                        #only be saved if the 'include_land' parameter in the
                        #sampling subdict is True);
                        #currently valid values: 'geotiff', 'txt'
                },
            }, #<END> 'data'

        'stats': {
            #dictionary defining which stats to be calculated, and parameters for 
                #their calculation (including frequency, in timesteps, of collection)
                #valid stats include:
                    # 'Nt'  : population census size
                    # 'het' : heterozygosity
                    # 'maf' : minor allele frequency
                    # 'ld'  : linkage disequilibrium
            'Nt':       {'calc' : True, 
                         'freq': 1,
                        },
            'het':      {'calc' : True, 
                         'freq': 1,
                         'mean': False,
                        },
            'maf':      {'calc' : True, 
                         'freq': 1,
                        },
            'ld':       {'calc' : True, 
                         'freq': 1,
                        }
            }, # <END> 'stats'

        } # <END> 'model'

    } # <END> 'params'



#a dict class with k:v pairs as dynamic attributes
class _Dyn_Attr_Dict_(dict):
    def __getattr__(self, item):
        return self[item]
    def __dir__(self):
        return super().__dir__() + [str(k) for k in self.keys()]
    def __deepcopy__(self, memo):
        return _Dyn_Attr_Dict_(copy.deepcopy(dict(self)))


class Params_Dict(_Dyn_Attr_Dict_):
    def __init__(self, params):
        params_dict = make_params_dict(params)
        self.update(params)
    #re-enable deepcopy, because the class inherits from a dict
    def __deepcopy__(self, memo):
        return Params_Dict(copy.deepcopy(dict(self)))


#function to recurse over the params dictionary 
#and return it as a Params_Dict object (i.e. a
#dict with k:v pairs as dynamic attributes)
def make_params_dict(params):
    for k, v in params.items():
        method_names = ['clear', 'copy', 'fromkeys', 'get', 'items', 'keys', 'pop', 'popitem', 'setdefault', 'update', 'values']
        assert k not in method_names, 'ERROR: The key "%s" in your params file is disallowed because it would clobber a Python method. Please edit name.\n\tNOTE: It holds the following value:\n%s' % (str(k), str(v))
        if isinstance(v, dict):
            #params.update({k:params_dict(v)})
            params[k] = make_params_dict(params[k])
    params = _Dyn_Attr_Dict_(params)
    return(params)
   


#get the params as a Params_Dict object
params = Params_Dict(params)





# <END> params.py
