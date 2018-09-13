#!/usr/bin/python
#genome.py


'''
##########################################

Module name:              genome

Module contents:          - definition of the Trait, RecombinationPaths, and  GenomicArchitecture classes
                          - function for simulation of genomic architecture, simulation of a new genome, and associated functions


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               12-28-15
Documentation:            URL


##########################################
'''

#geonomics imports
from ops import mutation

#other imports
import numpy as np
import pandas as pd
from numpy import random as r
from collections import OrderedDict as OD
import random
import bitarray


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

class Trait:
    def __init__(self, idx, name, phi, n_loci, mu, scape_num, mean_alpha_dist, std_alpha_dist, gamma, univ_advant):
        self.idx = idx
        self.name = name
        self.phi = phi
        self.n_loci = n_loci
        if mu is None:
            mu = 0
        self.mu = mu
        self.scape_num = scape_num
        self.mean_alpha_dist = mean_alpha_dist
        self.std_alpha_dist = std_alpha_dist
        self.gamma = gamma
        self.univ_advant = univ_advant

        self.loci = np.int64([])
        self.alpha = np.array([])


    def get_phi(self, pop):
        if type(self.phi) in (float, int):
            phi = np.array([self.phi]*len(pop))
        else:
            phi = self.phi[pop.cells[:,1], pop.cells[:,0]]
        return(phi)

    def set_loci(self, loci):
        self.loci = np.hstack((self.loci, np.array([*loci])))
        self.n_loci = self.loci.size

class RecombinationPaths:
    def __init__(self, recomb_paths):
        self.recomb_paths = recomb_paths

    def get_paths(self, n):
        return(random.sample(self.recomb_paths, n))


class GenomicArchitecture:
    def __init__(self, p, h, r, g_params):
        self.x = 2              #ploidy (NOTE: for now will be 2 by default; later could consider enabling polyploidy)
        self.L = g_params.L              #total length (i.e. number of markers)
        self.l_c = g_params.l_c          #length of each chromosome
        self.p = p              #Dict of allele frequencies for the 1-alleles for all (numbered by keys) chromosomes in haploid genome
        self.pleiotropy  = g_params.pleiotropy  #True/False regarding whether to allow a locus to affect the phenotype of more than one trait; defaults to False
        self.h = h          #Dict of heterozygous effects for all loci, for all chroms
        self.sex = g_params.sex
        self.r = r              #Dict of recombination rates between each locus and the next, for all chroms (NOTE: first will be forced to 1/float(x), to effect independent segregation of chroms, after which recomb occurs as a crossing-over path down the chroms
        self.recomb_lookup_array_size = g_params.recomb_lookup_array_size
        self.n_recomb_paths = g_params.n_recomb_paths

        self.recomb_paths = None  #The recombination-paths object will be assigned here; used to speed up large
                                #quantities of binomial draws needed for recombination

        #attributes for neutral loci
        self.mu_neut = g_params.mu_neut            #genome-wide neutral mutation rate  
        self.neut_loci = set(range(self.L))     #a set to keep track of all loci that don't influence the phenotypes of any trait; defaults to all loci, then will be updated
        self.nonneut_loci = set() #a set to keep track of all loci that influence the phenotype of at least one 
                                      #trait; after burn-in, will be updated

        #attributes for deleterious loci
        self.mu_delet = g_params.mu_delet                #genome-wide deleterious mutation rate
        self.delet_loci = OD()
        self.shape_delet_s_dist = g_params.shape_delet_s_dist
        self.scale_delet_s_dist = g_params.scale_delet_s_dist

        #attribute for trait (i.e. adaptive) loci
        self.traits = None
        print([*g_params])
        if 'traits' in [*g_params]:
            self.traits = make_traits(g_params.traits)

        #other attributes associated with mutation
            #a set to contain all mutable loci (i.e. loci where mutations can land)
        self.mutable_loci = set()    #a set containing eligible loci for mutation; after burn-in, will be updated
            #set self.mu_tot, the total per-site, per-generation mutation rate
        mus = [mu for mu in (self.mu_neut, self.mu_delet) if mu is not None]
        if self.traits is not None:
            mus = mus + [trt.mu for trt in self.traits.values()]
        self.mu_tot = sum(mus)
        self.mut_fns = self.make_mut_fns_dict()

    #method to make a mut_fns dict, containing a function for each type of mutation for this population
    def make_mut_fns_dict(self):
        mut_fns = {}
        if self.mu_neut > 0:
            fn = lambda pop,offspring: mutation.do_neutral_mutation(pop, offspring)
            mut_fns.update({'neut': fn})
        if self.mu_delet > 0:
            fn = lambda pop,offspring: mutation.do_deleterious_mutation(pop, offspring)
            mut_fns.update({'delet': fn})
        if self.traits is not None:
            for trait_num in self.traits:
                if self.traits[trait_num].mu > 0:
                    fn = lambda pop,offspring: mutation.do_trait_mutation(pop, offspring, trait_num)
                    mut_fns.update({'t%i' % trait_num: fn})
        return mut_fns

    #method to draw mutation types for any number of mutations chosen to occur in a given timestep 
    def draw_mut_types(self, num):
        type_dict = {'neut': self.mu_neut,
                 'delet': self.mu_delet,
                **{'t%i' % (k):v.mu for k,v in self.traits.items()} }
        types = []
        probs = []
        for k,v in type_dict.items():
            types.append(k)
            probs.append(v)
        probs = [p/sum(probs) for p in probs]
        choices = r.choice(types, p = probs, size = num, replace = True)
        return(choices)

    #method for drawing an effect size for one or many loci 
    def draw_trait_alpha(self, trait_num, n=1):
        alpha = r.normal(self.traits[trait_num].mean_alpha_dist, self.traits[trait_num].std_alpha_dist, n)
        #set all effects to positive if the trait is monogenic (because effects will be added to 0)
        if self.traits[trait_num].n_loci == 1:
            alpha = np.abs(alpha)
        return(alpha)

    #method for drawing new trait or deleterious loci (for a mutation) from the currently neutral loci
    def draw_mut_loci(self):
        loci = r.choice([*self.mutable_loci])
        #TODO: SHOULD I INCLUDE HERE A STATEMENT THAT CHECKS IF len(gen_arch.mutable_loci) IS <= SOME SMALL
        #VALUE, AND IF SO FINDS FIXED SITES IN THE POPULATION AND ADDS THEM TO IT??? (so that I'm not running
        #that time-intensive procedure often, but to make sure I'm almost certain not to run out of mutable
        #sites (unless I by chance draw more mutations in the next turn than there are remaining mutables))
        #DOES THIS IMPLY ANY STRANGE, POTENTIALLY PROBLEMATIC POPGEN ARTEFACTS?
        return(loci)

    #method for drawing new deleterious mutational fitness effects
    def draw_delet_s(self):
        s = r.gamma(self.shape_delet_s_dist, self.scale_delet_s_dist)
        s = min(s, 1)
        return(s)

    #a method to set the heterozygosity values at certain loci
    def set_heterozgosity(self, loci, het_value):
        pass

    #method for assigning loci to traits 
    def set_trait_loci(self, trait_num, mutational=False, loci=None, alpha=None):
        #if this is not the result of a point mutation, but instead either an initial setup or manually
        #introduced, then grab the number of loci to be assigned
        if mutational == False:
            n = self.traits[trait_num].n_loci
        #otherwise, assign a single locus
        else:
            n = 1
        if loci is not None:
            if not np.iterable(loci):
                loci = [loci]
            loci = set([*loci])
        #else, draw loci randomly, either allowing pleiotropy or not
        elif not self.pleiotropy:
            loci = set(r.choice([*self.neut_loci], size = n, replace = False))
        elif self.pleiotropy:
            loci = set(r.choice(range(self.L), size = n, replace = False))
        #update the trait's loci
        self.traits[trait_num].set_loci(loci)
        #add these loci to self.non-neutral and remove from self.neut_loci, to keep track of all loci underlying any traits (for purpose of avoiding pleiotropy)
        self.nonneut_loci.update(loci)
        self.neut_loci.difference_update(loci)
        #if the effect size(s) is/are provided, use those
        if alpha is not None:
            if not np.iterable(alpha):
                alpha = np.array([alpha]) 
            effects = np.array([*alpha])
        #else, draw effects from a Gaussian dist with mean 0 and sigma provided by trait params (as per e.g. Yeaman and Whitlock 2011)
        else:
            effects = self.draw_trait_alpha(trait_num, n)
        #check that loci and effects are of equal length
        assert len(loci) == len(effects), 'Lengths of the two arrays containing the new trait loci and their effects are not equal.'
        #then add the effects to the trait's alpha array
        self.traits[trait_num].alpha = np.hstack((self.traits[trait_num].alpha, effects))

    
    #method for creating and assigning the r_lookup attribute
    def make_recomb_paths(self):
        self.recomb_paths = RecombinationPaths(make_recomb_paths_bitarrays(self))


    #method for showing all allele frequencies for the population
    def show_freqs(self, pop):
        populome = np.hstack([ind.genome for ind in pop.values()])
        freqs = populome.sum(axis = 1)/(2*populome.shape[0])
        plt.plot(range(self.L), self.p, ':r')
        plt.plot(range(self.L), freqs, '_b')
        plt.show()


    #method for pickling a genomic architecture
    def write_pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

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


def make_traits(traits_params):
    params_copy = {**traits_params}
    #pop out the number of traits to create
    num_traits = len(params_copy)
    #then for each of i traits, unpack the ith components of the remaining params to create the trait dict
    traits = {k: Trait(k, **v) for k,v in params_copy.items()}
    return(traits)


#simulate linkage values
def draw_r(g_params, recomb_rate_fn = None):
    
    #use custom recomb_fn, if provided
    if recomb_rate_fn != None:
        recomb_array = np.array([max(0,min(0.5, recomb_rate_fn())) for locus in range(L)])
        return(recomb_array)

    #otherwise, use default function with either default or custom param vals
    else: 
        L = g_params.L

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


#carry out recombination, using the lookup array in a GenomicArchitecture object
def make_recombinants(r_lookup, n_recombinants):
    recombinants = np.array([r.choice(r_lookup[i,], size = n_recombinants, replace = True) for i in range(len(r_lookup))])
    recombinants = np.cumsum(recombinants, axis = 0)%2
    return(recombinants)


def make_recomb_array(g_params, recomb_values):
    #get L (num of loci) and l_c (if provided; num of loci per chromsome) from the genome params dict
    L = g_params.L
    if ('l_c' in g_params.keys() and g_params['l_c'] != None and len(g_params['l_c']) > 1):
        l_c = g_params.l_c
        #and if l_c provided, check chrom lenghts sum to total number of loci
        assert sum(l_c) == L, 'The chromosome lengths provided do not sum to the number of loci provided.'
    else:
        l_c = [L]

    #if g_params.recomb_array (i.e a linkage map) manually provided (will break if not a list, tuple, or np.array), 
    #then set that as the recomb_array, and check that len(recomb_array) == L
    if recomb_values is not None:
        assert len(recomb_values) == L, ('The length of the the table is the '
        'custom genomic-architecture file must be equal to the stipulated '
        "genome length ('L').")
        recomb_array = recomb_values[:]

    #otherwise, create recomb array
    else:
        #if a custom recomb_fn is provided, grab it
        if 'recomb_rate_custom_fn' in g_params.values():
            if g_params['recomb_rate_custom_fn'] is not None:
                recomb_rate_fn = g_params.recomb_rate_custom_fn
                assert callable(recomb_rate_fn), "The 'recomb_rate_custom_fn' provided in the parameters appears not to be defined properly as a callable function."
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
    #NOTE: #Always set the first locus r = 0.5, to ensure independent assortment of homologous chromosomes
    recomb_array[0] = 0.5

    return(recomb_array, sorted(l_c))


#function to create a lookup array, for raster recombination of larger numbers of loci on the fly
    #NOTE: size argument ultimately determines the minimum distance between probabilities (i.e. recombination rates)
    #that can be modeled this way
def make_recomb_paths_bitarrays(genomic_architecture, lookup_array_size = 10000, n_recomb_paths = 100000):
    
    if genomic_architecture.recomb_lookup_array_size is not None:
        lookup_array_size = genomic_architecture.recomb_lookup_array_size

    if genomic_architecture.n_recomb_paths is not None:
        n_recomb_paths = genomic_architecture.n_recomb_paths
    
    lookup_array = np.zeros((len(genomic_architecture.r),lookup_array_size), dtype = np.int8)

    for i, rate in enumerate(genomic_architecture.r):
        lookup_array[i,0:int(round(lookup_array_size*rate))] = 1

    recomb_paths = make_recombinants(lookup_array, n_recomb_paths).T
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
def make_genomic_architecture(pop_params):

    #get the genome parameters
    g_params = pop_params.gen_arch

    #get the custom genomic-architecture file, if provided
    gen_arch_file = None
    if 'gen_arch_file' in g_params.keys():
        if g_params.gen_arch_file is not None:
            gen_arch_file = pd.read_csv(g_params.gen_arch_file)
            assert len(gen_arch_file) == g_params.L, ('The custom '
            'genomic architecture file must contain a table with number of '
            'rows equal to the genome length stipulated in the genome '
            'parameters (params.comm[<pop_num>].gen_arch.L).')
            assert np.all(gen_arch_file['locus'].values == 
                          np.array([*range(len(gen_arch_file))])), ('The '
            "'locus' column of the custom genomic architecture file must "
            "contain serial integers from 0 to 1 - the length of the table.")

    #also get the sex parameter and add it as an item in g_params
    g_params['sex'] = pop_params.mating.sex

    #draw locus-wise 1-allele frequencies, unless provided in custom gen-arch file
    if gen_arch_file is None: 
        p = draw_allele_freqs(g_params.L)
    else:
        p = gen_arch_file['p'].values

    #get the custom recomb_values, if provided in a custom gen-arch file
    recomb_values = None
    if gen_arch_file is not None:
        recomb_values = gen_arch_file['r'].values
        assert recomb_values[0] == 0.5, ("The top recombination rate value "
        "in the custom genomic-architecture file's 'r' column must be set "
        "to 0.5, to effect independent assortment of sister chromatids.")

    #draw locus-wise heterozygosity-effect values
        #TODO: THIS WHOLE FUNCTIONALITY NEEDS TO BE EITHER RETOOLED OR ELSE ERADICATED
    h = draw_h(g_params.L)
    
    r, l_c = make_recomb_array(g_params, recomb_values = recomb_values)
    #in case g_params.l_c was missing or None, because only a single chromosome is being simulated, then
    #replace g_params.l_c with the returned l_c value
    g_params.l_c = l_c

    #now make the gen_arch object
    gen_arch = GenomicArchitecture(p,h, r, g_params)

    #set the loci and effect sizes for each trait, using the custom gen-arch
    #file, if provided
    if gen_arch_file is not None: 
        trait_effects = {trait_num: 
            {int(k): v for k,v, in gen_arch_file[gen_arch_file['trait'] == trait_num][['locus','alpha']].values}
            for trait_num in gen_arch.traits.keys()}
        #add the loci and effect sizes for each of the traits
        for trait_num in gen_arch.traits.keys():
            gen_arch.set_trait_loci(trait_num, mutational = False, 
                                loci = trait_effects[trait_num].keys(), 
                                alpha = trait_effects[trait_num].values())
    #or else randomly set the loci and effect sizes for each trait
    else:
        if gen_arch.traits is not None:
            for trait_num in gen_arch.traits.keys():
                gen_arch.set_trait_loci(trait_num, mutational = False)

    assert len(set(range(gen_arch.L)).difference(gen_arch.neut_loci.union(gen_arch.nonneut_loci))) == 0, 'ERROR: The union of the gen_arch.neut_loci and gen_arch.nonneut_loci sets does not contain all loci indicated by gen_arch.L'

    #create the r_lookup attribute
    gen_arch.make_recomb_paths()

    return(gen_arch)


#make genome
def draw_genome(genomic_architecture):
    new_genome = np.ones([genomic_architecture.L, genomic_architecture.x], dtype = np.int8)*9 #if for some reason any loci are not properly set to either 0 or 1, they will stick out as 9s
    for homologue in range(genomic_architecture.x):
        new_genome[:,homologue] = draw_genotype(genomic_architecture.p)

    assert type(new_genome) == np.ndarray, "A new genome must be an instance of numpy.ndarray"
    assert np.shape(new_genome) == (genomic_architecture.L, genomic_architecture.x), "A new genome must wind up with shape = (L, ploidy)."

    return(new_genome)


#function to reset genomes after burn-in
def set_genomes(pop, burn_T, T):
    from ops import mutation

    #use mean n_births at tail end of burn-in to estimate number of mutations, and randomly choose set of neutral loci 
    #of that length to go into the pop.gen_arch.mutable_loci attribute
    n_muts = mutation.calc_estimated_total_mutations(pop, burn_T, T)
    #add a small number if n_muts evaluates to 0
    if n_muts == 0:
        neut_loci_remaining = pop.gen_arch.L - len(pop.gen_arch.nonneut_loci)
        n_muts = min(int(np.ceil(0.1*neut_loci_remaining)), neut_loci_remaining)
    muts = set(r.choice([*pop.gen_arch.neut_loci], n_muts, replace = False))
    pop.gen_arch.mutable_loci.update(muts)

    #set those loci's p values to 0 (i.e. non-segregating)
    pop.gen_arch.p[np.array([*muts])] = 0
    
    #now reassign genotypes to all individuals, using gen_arch.p
    [ind.set_genome(draw_genome(pop.gen_arch)) for ind in pop.values()]
    #and then reset the individuals' phenotypes
    if pop.gen_arch.traits is not None:
        [ind.set_phenotype(pop.gen_arch) for ind in pop.values()];

       
#method for loading a pickled genomic architecture
def read_pickled_genomic_architecture(filename):
    import cPickle
    with open(filename, 'rb') as f:
        gen_arch = cPickle.load(f)

    return gen_arch

