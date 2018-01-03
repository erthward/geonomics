#!/usr/bin/python
#genome.py


'''
##########################################

Module name:              genome

Module contents:          - definition of the Genome and Genomic_Architecture classes
                          - function for simulation of genomic architecture, simulation of a new genome, and associated functions


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               12-28-15
Documentation:            URL


##########################################
'''


import numpy as np    
from numpy import random as r




#------------------------------------
# CLASSES ---------------------------
#------------------------------------


    

class Genomic_Architecture:
    def __init__(self, trait_dict, x, n, L, l_c, p, s, h, pleiotropy, env_var, r, mu, sex):
        self.x = x              #ploidy (NOTE: for now will be 2 by default; later could consider enabling polyploidy)
        self.n = n              #haploid number of chromosomes
        self.L = L              #total length (i.e. number of markers)
        self.l_c = l_c          #length of each chromosome
        self.p = p              #Dict of allele frequencies for the 1-alleles for all (numbered by keys) chromosomes in haploid genome
        self.s = s              #Dict of selection coefficients, for all chroms
        self.pleiotropy  = pleiotropy  #True/False regarding whether to allow a locus to affect the phenotype of more than one trait; defaults to False
        self.non_neutral = np.zeros([L]) #array to keep track of all loci that influence the phenotype of at least one trait
        self.h = h          #Dict of heterozygous effects for all loci, for all chroms
        self.env_var = env_var  #Dict of the environmental variables for which all loci on all chroms are selective
        self.r = r              #Dict of recombination rates between each locus and the next, for all chroms (NOTE: first will be forced to 1/float(x), to effect independent segregation of chroms, after which recomb occurs as a crossing-over path down the chroms
        self.mu = mu            #genome-wide mutation rate  #NOTE: should I allow this to be declared as a dictionary of mutation rates along all the chromosomes, to allow for heterogeneous mutation rates across the genome???
        self.sex = True 
        self.traits = trait_dict

    
   
    #method for drawing an effect size to a locus (or optionally, multiple)
    def draw_alpha(self, trait_num, n = 1):
        return (r.normal(0, self.traits[trait_num]['alpha_dist_sigma'], n))

    
    #method for assigning loci to traits 
    def assign_loci_to_trait(self, trait_num, mutational = False, locs = None):

        if mutational == False:  #i.e. if this is not the result of a point mutation, but instead either an initial setup or somehow manually introduced
            #number of loci to be assigned
            n = self.traits[trait_num]['num_loci']
        else:
            n = 1

        #if locations provided manually, use those
        if locs <> None:
            loci = np.array(locs)
            assert len(np.shape(loci)) == 1, 'Object containing loci provided manually is not of dimension 1'
        #otherwise choose loci randomly, either allowing pleiotropy or not
        elif self.pleiotropy == False:
            loci = r.choice(np.where(self.non_neutral == False)[0], size = n, replace = False)
        elif self.pleiotropy == True:
            loci = r.choice(range(self.L), size = n, replace = False)

        #choose effects from a Gaussian dist with mean 0 and sigma provided by trait params (as per e.g.  Yeaman and Whitlock 2011)
        effects = self.draw_alpha(trait_num, n)

        #check that loci and effects are of equal length
        assert len(loci) == len(effects), 'Lengths of the two arrays containing the new trait loci and their effects are not equal.'

        #add these loci to the trait's 'loci' array
        self.traits[trait_num]['loci'] = np.array(list(self.traits[trait_num]['loci']) + list(loci))
        #and add their effects to the 'alpha' array
        self.traits[trait_num]['alpha'] = np.array(list(self.traits[trait_num]['alpha']) + list(effects))

        #and and these loci to self.non-neutral, to keep track of all loci underlying any traits (for purpose of avoiding pleiotropy)
        self.non_neutral[loci] = 1



    #method for pickling a genomic architecture
    def pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)





class Genome:
    def __init__(self, G, ploidy):
        self.genome = G #dict of numbered i*j arrays, each subarray containing the genotype data for all j (for a diploid, 2) copies of chromosome i
        assert list(set([type(a) for a in self.genome.values()])) == [type(np.array([0]))], "G must be a dict of numpy.ndarray instances"
        assert list(set([np.shape(a)[1] for a in self.genome.values()])) == [ploidy] 

    #NOTE: need to ensure that size(G)[1] == n and that size of each subarray in G == x*l_c




#----------------------------------
# FUNCTIONS -----------------------
#----------------------------------


#function for choosing a random number from a Bernoulli distribution of an arbitrary p value
#(to be mapped to allele frequencies in order to generate genotypes)
def rand_bern(prob):
    return r.binomial(1,prob)



#generate allele_freqs
def gen_allele_freqs(l):
    return r.beta(1,1,l)



#simulate genotypes
def sim_G(p): 
    return [rand_bern(freq) for freq in p]




#randomly assign heterozygous effects, choosing from either 0, 0.5, or 1,
#where 0 = allele 1 (A1) dominant, 0.5 = codominant, 1 = A1 recessive, i.e.
#relative fitnesses: A1A1 = 1, A1A2 = 1-hs, A2A2 = 1-s
#NOTE: Should be able to also implement negative h values here, for overdominance
def assign_h(l, params, low = 0, high = 2):
    h  = r.randint(low = 0, high = 3, size = l)/2.   #NOTE: Pythonic style, so high is exclusive, hence high = 3
    return h




#randomly assign selectivity of each locus to one of the eligible rasters 
    #NOTE: Original idea was that 0 = global selection (i.e.  selective advantage non-variable in space, i.e.
    #fitness of alleles will not be contingent on any raster), but got rid of that on 11/19/17
def assign_env_var(num_rasters, l, params, allow_multiple_env_vars = True):


    #if allow_multiple_env_vars, then highest env_var number will be the last landscape in the landscape_stack 
    if allow_multiple_env_vars == True:
        high = num_rasters
    else:
        high = 1

    possible_vars = range(high)

    if ('movement_surf' in params.keys()):
        if params['movement_surf'] == True:
            exclude_var = params['movement_surf_scape_num']
            possible_vars.remove(exclude_var)

    
    env_var = r.choice(possible_vars, size = l, replace = True)


    return env_var




  
#simulate selection coefficients
def sim_s(l, alpha_s = 0.007, beta_s = 2): #NOTE: alpha= 0.15 gives ~600-660 loci in 10,000 with s > .75; 0.025 gives ~85-115; 0.0025 gives ~5-15
    #NOTE: See Thurman and Barrett (2016) for metaanalysis of s values in real populations!

    s = r.beta(alpha_s, beta_s, l)
    s = np.array([0.0 if l < 0.001 else l for l in s])
    return s
    #NOTE: For now, beta seems an intuitive and fine way to model this
    #For effectively no selection: Beta(1,1e5)
    #For spread of different selection values, but highly zero-inflated: ~ Beta(0.15, 2)
    #For a wide range of selection coefficients, something like: Beta(1,5)
    #For a uniform range of selection coefficients between 0 and 1: beta (1,1)




def construct_trait_dict(traits_params):
    trait_dict = {}
    #for each trait, create a dict with the params
    for t in range(traits_params['num']):
        curr_trait_dict = dict([(k,v[t]) for k,v in traits_params.items() if k <> 'num'])
        #use the params to create additional necessary dict components, including:
        #loci assigned to the trait (initially empty, will be filled by calling a Genomic_Architecture method after genomic_arch is created)
        curr_trait_dict['loci'] = []
        #effect sizes of those loci (also initally empty)
        curr_trait_dict['alpha'] = []
        trait_dict[t] = curr_trait_dict
    return(trait_dict)



#simulate linkage values
def sim_r(l, alpha_r = 7e2, beta_r = 7e3):
    return np.array([min(0.5, n) for n in r.beta(alpha_r, beta_r, l)])
    #NOTE: for now, using the Beta, which can be very flexibly parameterized
    #NOTE: current default alpha/beta vals, after being subtracted from 0.5 in sim_r function, will result in a tight distribution of r vals around 0.21 (~ 95% between 0.19 and 0.22)
    #NOTE: to essentially fix r at 0.5, Beta(1,1e7) will do...
    #NOTE: Ideally, would be good to look into to developing some sort of mixture distribution to reasonably and flexibly model map-distance values...



#set chromosome lengths
def set_l_c(L, n, even_chrom_sizes = True):
    if not even_chrom_sizes:
        None #NOTE: Instead, if command-line arg provided, allow it to use user-supplied array stipulating of chromosome lengths
        return l_c
    else:
        l_c = [int(round(L/n)) for i in range(n)]
        return l_c



#build the genomic architecture
#NOTE: This will create the "template" for the genomic architecture that will then be used to simulate individuals and populations
def build_genomic_arch(params, land, allow_multiple_env_vars = True):




    #NOTE: 11/19/17: DO I STILL WANT TO OPERATIONALIZE GLOBALLY SELECTIVE LOCI??


    #grab necessary parameters from the params dict

    L = params['L']
    n = params['n']
    mu = params['mu']
    x = params['x']
    alpha_r = params['alpha_r']
    beta_r = params['beta_r']
    pleiotropy = params['pleiotropy']

    trait_dict = construct_trait_dict(params['traits'])




    #NOTE: THIS SEEMS LIKE A VESTIGE FROM SOME PREVIOUS IDEA THAT IS NOW NOT CLEAR TO ME... INVESTIGATE, THEN LIKELY TEAR OUT
    sex = params['sex']  #NOTE: HOW TO CHANGE THIS TO MAKE USE OF THE 'sex' PARAM IN params???
    if True:
        sex = False



    l_c = None   #NOTE: OPERATIONALIZE UNEVEN CHROM LENGTHS??





    if l_c:
        pass #Allow provision of array of uneven chromosome lengths
    else:
    #NOTE: x = ploidy, for now set to 2 (i.e.  diploidy)
    #NOTE: how to operationalize sexuality?! for now defaults to False
        l_c = set_l_c(L, n)
        p = dict()
        s = dict()
        h = dict()
        env_var = dict()
        r = dict()
        for chrom in range(n):
            p[chrom] = gen_allele_freqs(l_c[chrom])
            s[chrom] = sim_s(l_c[chrom])
            h[chrom] = assign_h(l_c[chrom], params)
            env_var[chrom] = assign_env_var(num_rasters = land.num_rasters, l = l_c[chrom], params = params, allow_multiple_env_vars = allow_multiple_env_vars)
            r[chrom] = sim_r(l_c[chrom], alpha_r, beta_r)
            r[chrom][0] = 1/float(x) #setting first val to 1/float(x) effects independent segregation of chroms, after which recombination occurs along the chrom according to defined recomb rates (r)
        genomic_arch = Genomic_Architecture(trait_dict, x, n, sum(l_c), l_c, p, s, h, pleiotropy, env_var, r, mu, False)
        
        #add the loci and effect sizes for each of the traits
        for trait_num in genomic_arch.traits.keys():
            genomic_arch.assign_loci_to_trait(trait_num, mutational = False, locs = None)

        return Genomic_Architecture(trait_dict, x, n, sum(l_c), l_c, p, s, h, pleiotropy, env_var, r, mu, False)





#simulate genome
def sim_genome(genomic_arch):
    new_genome = dict()
    for chrom in range(genomic_arch.n):
        chromosome = np.ones([genomic_arch.l_c[chrom], genomic_arch.x])*999 #if for some reason any loci are not properly set to either 0 or 1, they will stick out as 999's
        for copy in range(genomic_arch.x):
            chromosome[:,copy] = sim_G(genomic_arch.p[chrom])
        new_genome[chrom] = chromosome
    return Genome(new_genome, genomic_arch.x)




#method for loading a pickled genomic architecture
def load_pickled_genomic_arch(filename):
    import cPickle
    with open(filename, 'rb') as f:
        genomic_arch = cPickle.load(f)

    return genomic_arch


