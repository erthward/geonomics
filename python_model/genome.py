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
import random

import bitarray
import gametogenesis




#------------------------------------
# CLASSES ---------------------------
#------------------------------------


class Trait:
    def __init__(self, num, phi, n_loci, scape_num, alpha_dist_sigma, gamma, univ_advant):
        self.num = num
        self.phi = phi
        self.n_loci = n_loci
        self.scape_num = scape_num
        self.alpha_dist_sigma = alpha_dist_sigma
        self.gamma = gamma
        self.univ_advant = univ_advant

        self.loci = np.array([])
        self.alpha = np.array([])


    def get_phi(self, pop):
        if type(self.phi) in (float, int):
            phi = np.array([self.phi]*pop.census())
        else:
            phi = self.phi[pop.cells[:,1], pop.cells[:,0]]
        return(phi)


class Recomb_Paths:
    def __init__(self, recomb_paths):
        self.recomb_paths = recomb_paths

    def get(self, n):
        return(random.sample(self.recomb_paths, n))


class Genomic_Architecture:
    def __init__(self, x, L, l_c, r, mu, p, trait_dict, h, pleiotropy, sex):
        self.x = x              #ploidy (NOTE: for now will be 2 by default; later could consider enabling polyploidy)
        self.L = L              #total length (i.e. number of markers)
        self.l_c = l_c          #length of each chromosome
        self.p = p              #Dict of allele frequencies for the 1-alleles for all (numbered by keys) chromosomes in haploid genome
        self.pleiotropy  = pleiotropy  #True/False regarding whether to allow a locus to affect the phenotype of more than one trait; defaults to False
        self.non_neutral = np.zeros([L]) #array to keep track of all loci that influence the phenotype of at least one trait
        self.h = h          #Dict of heterozygous effects for all loci, for all chroms
        self.r = r              #Dict of recombination rates between each locus and the next, for all chroms (NOTE: first will be forced to 1/float(x), to effect independent segregation of chroms, after which recomb occurs as a crossing-over path down the chroms
        self.recomb_paths = None  #The recombination-paths object will be assigned here; used to speed up large
                                #quantities of binomial draws needed for recombination
        self.mu = mu            #genome-wide mutation rate  #NOTE: should I allow this to be declared as a dictionary of mutation rates along all the chromosomes, to allow for heterogeneous mutation rates across the genome???
        self.mutables = None    #after burn-in, will be set to an array containing eligible loci for mutation
        self.sex = sex

        self.traits = trait_dict

    
   
    #method for drawing an effect size for one or many loci 
    def draw_alpha(self, trait_num, n):
        alpha = r.normal(0, self.traits[trait_num].alpha_dist_sigma, n)
        if n == 1:
            alpha = np.abs(alpha)
        return(alpha)

    
    #method for assigning loci to traits 
    def assign_loci_to_trait(self, trait_num, mutational = False, locs = None):

        if mutational == False:  #i.e. if this is not the result of a point mutation, but instead either an initial setup or somehow manually introduced
            #number of loci to be assigned
            n = self.traits[trait_num].n_loci
        else:
            n = 1

        #if locations provided manually, use those
        if locs != None:
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
        self.traits[trait_num].loci = np.array(list(self.traits[trait_num].loci) + list(loci))
        #and add their effects to the 'alpha' array
        self.traits[trait_num].alpha = np.array(list(self.traits[trait_num].alpha) + list(effects))

        #and add these loci to self.non-neutral, to keep track of all loci underlying any traits (for purpose of avoiding pleiotropy)
        self.non_neutral[loci] = 1


    #method for creating and assigning the r_lookup attribute
    def create_recomb_paths(self):
        self.recomb_paths = Recomb_Paths(create_recomb_paths_bitarrays(self))




    #method for pickling a genomic architecture
    def pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)







#----------------------------------
# FUNCTIONS -----------------------
#----------------------------------



#generate allele_freqs
def draw_allele_freqs(l):
    return(r.beta(1,1,l))



#simulate genotypes
def draw_genotype(p): 
    return(r.binomial(1, p))


#randomly assign heterozygous effects, choosing from either 0, 0.5, or 1,
#where 0 = allele 1 (A1) dominant, 0.5 = codominant, 1 = A1 recessive, i.e.
#relative fitnesses: A1A1 = 1, A1A2 = 1-hs, A2A2 = 1-s
#NOTE: Should be able to also implement negative h values here, for overdominance
def draw_h(l, low = 0, high = 2):
    h  = r.randint(low = 0, high = 3, size = l)/2.   #NOTE: Pythonic style, so high is exclusive, hence high = 3
    return(h)



def construct_traits(traits_params):
    from copy import deepcopy
    params_copy = deepcopy(traits_params)
    #pop out the number of traits to create
    num_traits = params_copy.pop('num') 
    #then for each of i traits, unpack the ith components of the remaining params to create the trait dict
    traits = {i : Trait(i, **{k:v[i] for k,v in params_copy.items()}) for i in range(num_traits)}
    return(traits)



#simulate linkage values
def draw_r(g_params, recomb_rate_fn = None):
    
    #use custom recomb_fn, if provided
    if recomb_rate_fn != None:
        recomb_array = np.array([max(0,min(0.5, recomb_rate_fn())) for locus in range(L)])
        return(recomb_array)

    #otherwise, use default function with either default or custom param vals
    else: 
        L = g_params['L']

        param_vals = {'alpha_r': 7e2, 'beta_r': 7e3}

        for param in ['alpha_r', 'beta_r']:
            if (param in g_params.keys() and g_params[param] != None):
                param_vals[param] = g_params[param]
       
        recomb_array = np.array([max(0,min(0.5, recomb_rate)) for recomb_rate in r.beta(a = param_vals['alpha_r'], b = param_vals['beta_r'], size = L)])
    #NOTE: for now, using the Beta, which can be very flexibly parameterized
    #NOTE: current default alpha/beta vals, after being subtracted from 0.5 in sim_r function, will result in a tight distribution of r vals around 0.21 (~ 95% between 0.19 and 0.22)
    #NOTE: to essentially fix r at 0.5, Beta(1,1e7) will do...
    #NOTE: Ideally, would be good to look into to developing some sort of mixture distribution to reasonably and flexibly model map-distance values...

        return(recomb_array)




def get_chrom_breakpoints(l_c, L):

    breakpoints = np.array([0]+list(np.cumsum(sorted(l_c))[:-1]))

    assert np.alltrue(np.diff(np.array(list(breakpoints) + [L])) == np.array(sorted(l_c))), 'The breakpoints assigned will not produce chromosomes of the correct length'

    return(breakpoints)




def create_recomb_array(g_params):


    #get L (num of loci) and l_c (if provided; num of loci per chromsome) from params['genome'] dict
    L = g_params['L']
    if ('l_c' in g_params.keys() and g_params['l_c'] != None and len(g_params['l_c']) > 1):
        l_c = g_params['l_c']
        #and if l_c provided, check chrom lenghts sum to total number of loci
        assert sum(l_c) == L, 'The chromosome lengths provided do not sum to the number of loci provided.'
    else:
        l_c = [L]


    #if g_params['recomb_array'] (i.e a linkage map) manually provided (will break if not a list, tuple, or np.array), 
    #then set that as the recomb_array, and check that len(recomb_array) == L
    if ('recomb_array' in g_params.keys() and g_params['recomb_array'] != None):
        recomb_array = np.array(g_params['recomb_array'])
        assert len(recomb_array) == L, "Length of recomb_array provided not equal to params['genome']['L']."

        #NOTE: #Always necessary to set the first locus r = 1/ploidy, to ensure independent assortment of homologous chromosomes
        recomb_array[0] = 0.5
        #NOTE: for now, obligate diploidy
        #recomb_array[0] = 1/g_params['x'] 

        return(recomb_array)

    
    #otherwise, create recomb array
    else:

        #if a custom recomb_fn is provided, grab it
        if ('recomb_rate_custom_fn' in g_params.values() and g_params['recomb_rate_custom_fn'] is not None):
            recomb_rate_fn = g_params['custom_fns']
            assert callable(recomb_rate_fn), "The recomb_rate_fn provided in params['genome']['recomb_rate_custom_fn'] appear not to be defined properly as a callable function."
            #then call the draw_r() function for each locus, using custom recomb_fn
            recomb_array = draw_r(g_params, recomb_fn = recomb_rate_fn)



        #otherwise, use the default draw_r function to draw recomb rates
        else:
            recomb_array = draw_r(g_params) 


        #if more than one chromosome (i.e. if l_c provided in g_params dict and of length >1), 
        #set recomb rate at the appropriate chrom breakpoints to 0.5
        if len(l_c) >1:
            bps = get_chrom_breakpoints(l_c, L)
            recomb_array[bps] = 0.5


        #NOTE: #Always necessary to set the first locus r = 0.5, to ensure independent assortment of homologous chromosomes
        recomb_array[0] = 0.5
        #NOTE: for now, obligate diploidy
        #recomb_array[0] = 1/g_params['x']


        return(recomb_array, sorted(l_c))




#function to create a lookup array, for raster recombination of larger numbers of loci on the fly
    #NOTE: size argument ultimately determines the minimum distance between probabilities (i.e. recombination rates)
    #that can be modeled this way
def create_recomb_paths_bitarrays(genomic_arch, lookup_array_size = 10000, num_recomb_paths = 10000):
    
    lookup_array = np.zeros((len(genomic_arch.r),lookup_array_size), dtype = np.int8)

    for i, rate in enumerate(genomic_arch.r):
        lookup_array[i,0:int(round(lookup_array_size*rate))] = 1

    recomb_paths = gametogenesis.recombine(lookup_array, num_recomb_paths).T
    bitarrays = tuple([make_bitarray_recomb_subsetter(p) for p in recomb_paths])
    
    return(bitarrays)




def make_bitarray_recomb_subsetter(recomb_path):
    ba = bitarray.bitarray(list(recomb_path.reshape((recomb_path.size,))))
    ba_inv = bitarray.bitarray(list(np.array(ba) == False))
    tot = []
    for i in range(len(ba)):
        tot.append(ba[i])
        tot.append(ba_inv[i])
    return(bitarray.bitarray(tot))



#build the genomic architecture
#NOTE: This will create the "template" for the genomic architecture that will then be used to simulate individuals and populations
def build_genomic_arch(g_params, land, allow_multiple_env_vars = True):

    g_params = params['genome']

    p_params = params['pop']


    #NOTE: 11/19/17: DO I STILL WANT TO OPERATIONALIZE GLOBALLY SELECTIVE LOCI??


    #grab necessary parameters from the g_params dict
    L = g_params['L']
    mu = g_params['mu']
    x = g_params['x']
    pleiotropy = g_params['pleiotropy']


    trait_dict = construct_traits(g_params['traits'])




    #NOTE: THIS SEEMS LIKE A VESTIGE FROM SOME PREVIOUS IDEA THAT IS NOW NOT CLEAR TO ME... INVESTIGATE, THEN LIKELY TEAR OUT
    sex = p_params['sex']  #NOTE: HOW TO CHANGE THIS TO MAKE USE OF THE 'sex' PARAM IN p_params???
    if True:
        sex = False



    #NOTE: x = ploidy, for now set to 2 (i.e.  diploidy)
    #NOTE: how to operationalize sexuality?! for now defaults to False

    p = draw_allele_freqs(L)  
        #TODO: DECIDE HOW TO OPERATIONALIZE MUTATIONS; PERHAPS WANT TO CREATE A BUNCH OF p=0 LOCI HERE, AS MUTATIONAL TARGETS FOR LATER?

    h = draw_h(L)

    r, l_c = create_recomb_array(g_params)

    genomic_arch = Genomic_Architecture(x, L, l_c, r, mu, p, trait_dict, h, pleiotropy, sex = True)
   

    #add the loci and effect sizes for each of the traits
    for trait_num in genomic_arch.traits.keys():
        genomic_arch.assign_loci_to_trait(trait_num, mutational = False, locs = None)

    #create the r_lookup attribute
    genomic_arch.create_recomb_paths()


    return(genomic_arch)





#simulate genome
def sim_genome(genomic_arch):
    new_genome = np.ones([genomic_arch.L, genomic_arch.x], dtype = np.int8)*9 #if for some reason any loci are not properly set to either 0 or 1, they will stick out as 9s
    for homologue in range(genomic_arch.x):
        new_genome[:,homologue] = draw_genotype(genomic_arch.p)

    assert type(new_genome) == np.ndarray, "A new genome must be an instance of numpy.ndarray"
    assert np.shape(new_genome) == (genomic_arch.L, genomic_arch.x), "A new genome must wind up with shape = (L, ploidy)."

    return(new_genome)




#function to reassign genomes after burn-in
def reassign_genomes(pop, params):
    import mutation

    #use mean n_births at tail end of burn-in to estimate number of mutations, and randomly choose set of neutral loci 
    #of that length to go into pop.genomic_arch.mutable slot
    n_muts = mutation.estimate_total_num_mutations(params, pop)
    muts = r.choice(np.where(pop.genomic_arch.non_neutral == 0)[0], n_muts, replace = False)
    pop.genomic_arch.mutables = list(muts)

    #set those loci's p values to 0 (i.e. non-segregating)
    pop.genomic_arch.p[muts] = 0
    
    #now reassign genotypes to all individuals, using genomic_arch.p
    for ind in pop.individs.keys():
        pop.individs[ind].genome = sim_genome(pop.genomic_arch)

       
#method for loading a pickled genomic architecture
def load_pickled_genomic_arch(filename):
    import cPickle
    with open(filename, 'rb') as f:
        genomic_arch = cPickle.load(f)

    return genomic_arch

