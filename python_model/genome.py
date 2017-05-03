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


#NOTE: shouldn't need to import these here, right? since they're imported in main.py...
import numpy as np    
from numpy import random as r




#------------------------------------
# CLASSES ---------------------------
#------------------------------------


    

class Genomic_Architecture:
    def __init__(self, x, n, L, l_c, p, s, h, env_var, D, mu, sex):
        self.x = x              #ploidy (NOTE: for now will be 2 by default; later could consider enabling polyploidy)
        self.n = n              #haploid number of chromosomes
        self.L = L              #total length (i.e. number of markers)
        self.l_c = l_c          #length of each chromosome
        self.p = p              #Dict of allele frequencies for the 1-alleles for all (numbered by keys) chromosomes in haploid genome
        self.s = s              #Dict of selection coefficients, for all chroms
        self.non_neutral = dict([(c, self.s[c] <> 0) for c in self.s.keys()])  #np.array for quick subsetting of only non-neutral loci, for use in selection and mutation modules
        self.sum_s =     sum([sum(self.s[c][self.non_neutral[c]]) for c in self.s.keys()])  #sum of all selection coefficients; to save on computation in weighted sum used in selection module
        self.h = h          #Dict of heterozygous effects for all loci, for all chroms
        self.env_var = env_var  #Dict of the environmental variables for which all loci on all chroms are selective
        self.D = D              #Dict of linkage disequilibria (i.e. map distances) between each locus and the next, for all chroms (NOTE: because of this definition, the last value is forced to 0)
        self.mu = mu            #genome-wide mutation rate  #NOTE: should I allow this to be declared as a dictionary of mutation rates along all the chromosomes, to allow for heterogeneous mutation rates across the genome???
        self.sex = True 



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




#randomly assign selectivity of each locus to one of the eligible rasters (note: 0 = global selection (i.e.  selective advantage non-variable in space, i.e. fitness of alleles will not be contingent on any raster)
def assign_env_var(num_rasters, l, allow_multiple_env_vars = False, allow_global_selection = True):


    #if allow_multiple_env_vars, then highest env_var number will be the last landscape in the landscape_stack; otherwise, all env_var values will be either 1 (for the first raster) or 0 (for global selection if allow_global_selection = True)
    if allow_multiple_env_vars == True:
        high = num_rasters
    else:
        high = 1


    env_var = r.random_integers(low = not allow_global_selection, high = high, size = l)  #NOTE: 0's will be assigned if allow_global_selection = True, else 1 will be lowest number)


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



#simulate linkage values
def sim_D(l, alpha_D = 7e2, beta_D = 7e3):
    return r.beta(alpha_D, beta_D, l) 
    #NOTE: for now, using the Beta, which can be very flexibly parameterized
    #NOTE: current default alpha/beta vals, after being subtracted from 0.5 in sim_D function, will result in a tight distribution of D vals around 0.21 (~ 95% between 0.19 and 0.22)
    #NOTE: to essentially fix D at 0.5, Beta(1,1e7) will do...
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
#NOTE: This will create the "template" for the genomic architecture of the hapolid genome that will then be used to simulate individuals and populations
def build_genomic_arch(params, land, allow_global_selection = False, allow_multiple_env_vars = False):




    #NOTE: OPERATIONALIZE GLOBALLY SELECTIVE LOCI??
    #NOTE: 04-09-17: Set allow_multiple_env_vars to default to False




    #grab necessary parameters from the params dict

    L = params['L']
    n = params['n']
    mu = params['mu']
    x = params['x']
    alpha_D = params['alpha_D']
    beta_D = params['beta_D']




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
        D = dict()
        for chrom in range(n):
            p[chrom] = gen_allele_freqs(l_c[chrom])
            s[chrom] = sim_s(l_c[chrom])
            h[chrom] = assign_h(l_c[chrom], params)
            env_var[chrom] = assign_env_var(num_rasters = land.num_rasters, l = l_c[chrom], allow_global_selection = allow_global_selection, allow_multiple_env_vars = allow_multiple_env_vars)
            D[chrom] = sim_D(l_c[chrom], alpha_D, beta_D)
            D[chrom][len(D[chrom])-1] = 0 #because each D value expresses linkage between that locus and the next, last value is forced to 0!
        return Genomic_Architecture(x, n, sum(l_c), l_c, p, s, h, env_var, D, mu, False)




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


