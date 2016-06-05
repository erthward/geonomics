#!/usr/bin/python 
#population.py

'''
##########################################

Module name:              population


Module contains:
                          - definition of the Population class
                          - function for creating a population of Individuals (incl. their genomes and associated data)
                          - associated functions


Author:                    Drew Ellison Hart
Email:                     drew.hart@berkeley.edu
Github:                    URL
Start date:                12-28-15
Documentation:             URL


##########################################
'''

import genome
import individual
import mating
import dispersal
import selection
import mutation

import numpy as np
from numpy import random as r
import matplotlib as mpl
import matplotlib.pyplot as plt

import sys

#------------------------------------
# CLASSES ---------------------------
#------------------------------------


class Population:



    '''Population(self, N, individs, genomic_arch, size, birth_rate = None, death_rate = None, T = None)'''


     
    def __init__(self, N, individs, genomic_arch, size, T): 

        self.N = []                                        #list to record population size at each step (with starting size as first entry)
        self.individs = individs                            #dict of all individuals, as instances of individual.Individual class


        self.initial_size = len(self.individs)


        
        if size.__class__.__name__ in ['float', 'int'] and T <> None:
            self.size = [size] * T              
        elif size.__class__.__name__ in ['list']:
            assert T<> None and len(size) == T, "if expressing population size as a list, total model runtime must be provided and must equal the list length"
            self.size = size
        else:
            self.size = size



        self.genomic_arch = genomic_arch
        #Add other demographic parameters?? Other stuff in general??



        assert type(N) == int, "N must be an integer"
        assert type(self.individs) == dict, "self.individs must be a dictionary"
        assert list(set([i.__class__.__name__ for i in self.individs.values()])) == ['Individual'], "self.individs must be a dictionary containing only instances of the individual.Individual class"
        assert self.genomic_arch.__class__.__name__ == 'Genomic_Architecture', "self.genomic_arch must be an instance of the genome.Genomic_Architecture class"





    #####################
    ### OTHER METHODS ###
    #####################

    def census(self):
        return len(self.individs)



    #method to increment all population's age by one (also adds current pop size to tracking array)
    def birthday(self):
        #add current pop size to pop.N (for later demographic analysis)
        self.N.append(self.census())
        #increment age of all individuals
        [ind.birthday() for ind in self.individs.values()];



    #method to extract habitat values at each individual's coordinates for each land.scapes raster
    def query_habitat(self, land):
        [ind.query_habitat(land) for ind in self.individs.values()];
 




    #method to move all individuals simultaneously
    def move(self, land, params, resist_surf = None):
        [ind.move(land, params['mu_direction'], params['kappa_direction'], params['mu_distance'], params['sigma_distance'], resist_surf) for ind in self.individs.values()];

        self.query_habitat(land)




    #set mating.find_mates() as method
    def court(self, params):

        return mating.find_mates(self, params['mating_radius'], sex = params['sex'], repro_age = params['repro_age'])





    #function for executing mating for a population
    def mate(self, land, params):

        
        
        #pull necessary parameters from params dict
        mating_radius = params['mating_radius']
        mu_dispersal = params['mu_dispersal']
        sigma_dispersal = params['sigma_dispersal']
        sex = params['sex']
        repro_age = params['repro_age']



        mating_pairs = self.court(params)

        num_offspring = 0

        for pair in mating_pairs:

            parent_centroid_x = np.mean((self.individs[pair[0]].x, self.individs[pair[1]].x))
            parent_centroid_y = np.mean((self.individs[pair[0]].y, self.individs[pair[1]].y))

            zygotes = mating.mate(self, pair, self.genomic_arch)

            
            for zygote in zygotes:

                num_offspring += 1
                
                offspring_x, offspring_y = dispersal.disperse(land, parent_centroid_x, parent_centroid_y, mu_dispersal, sigma_dispersal)

                if sex == True:
                    offspring_sex = r.binomial(1, 0.5)

                age = 0

                offspring_key = max(self.individs.keys()) + 1

                if sex == True:
                    self.individs[offspring_key] = individual.Individual(zygote, offspring_x, offspring_y, offspring_sex, age)
                else:
                    self.individs[offspring_key] = individual.Individual(zygote, offspring_x, offspring_y, age)


        #sample all individuals' habitat values, to initiate for offspring
        self.query_habitat(land)

        print '\t%i individuals born' % num_offspring






    #method to carry out selection
    def select(self, t, params):
        selection.select(self, t, params, sigma_deaths = params['sigma_deaths'], density_dependent_fitness = params['density_dependent_fitness'])






    #method to carry out mutation
    def mutate(self, params, t):
        for ind in [ind for ind, individ in self.individs.items() if individ.age == 0]:
            mutation.mutate(self, self.individs[ind], self.genomic_arch, t, alpha_mut_s = params['alpha_mut_s'], beta_mut_s = params['beta_mut_s'])





    def check_extinct(self):
        if len(self.individs.keys()) == 0:
            print '\n\nYOUR POPULATION WENT EXTINCT!\n\n\t(Press <Enter> to exit.)'
            raw_input()
            sys.exit()
    




    

    #function to discover loci above, below, or between threshold selection coefficients
    def find_loci(self, min_s=None, max_s = None):
        assert min_s<> None or max_s <> None, "No parameters provided. Must provide at least one value (either a max or min)."
        if max_s and not min_s:
            loci = dict([(i, np.array(range(len(self.genomic_arch.s[i])))[self.genomic_arch.s[i] <= max_s]) for i in range(len(self.genomic_arch.s))])
            #print "\n%i loci found\n" % sum([len(chrom_set) for chrom_set in loci.values()])
            return loci


        elif min_s and not max_s:
            loci = dict([(i, np.array(range(len(self.genomic_arch.s[i])))[self.genomic_arch.s[i] >= min_s]) for i in range(len(self.genomic_arch.s))])
            #print "\n%i loci found\n" % sum([len(chrom_set) for chrom_set in loci.values()])
            return loci
 

        else:
            loci = dict([(i, np.array(range(len(self.genomic_arch.s[i])))[np.array(self.genomic_arch.s[i] >= min_s) & np.array(self.genomic_arch.s[i] <= max_s)]) for i in range(len(self.genomic_arch.s))])
            #print "\n%i loci found\n" % sum([len(chrom_set) for chrom_set in loci.values()])
            return loci
 

   


   
 

    def get_habitat(self, individs = None):
        if individs <> None:
            return dict([(k, ind.habitat) for k, ind in self.individs.items() if k in individs])
        else:
            return dict([(k, ind.habitat) for k, ind in self.individs.items()])


    

    def get_genotype(self, chromosome, locus, format = 'mean', individs = None):

        dom_funcs = { 0 : np.mean,                          #codominance
                      1 : lambda x: np.ceil(np.mean(x)),    #dominant (with respect to relationship between allele 1 and habitat values -> 1
                      2 : lambda x: np.floor(np.mean(x))    #recessive (with respect to relationship between allele 1 and habitat values -> 1
                    }


        if individs == None:
            individs = self.individs.keys()
            #individs = range(len(self.genomic_arch.s[chromosome]))

        if format == 'biallelic':
            return dict([(i, self.individs[i].genome.genome[chromosome][locus, :]) for i in self.individs.keys() if i in individs]) 

        elif format == 'mean':
            return dict([(i, [dom_funcs[self.genomic_arch.dom[chromosome][locus]](self.individs[i].genome.genome[chromosome][locus,:]) ]) for i in self.individs.keys() if i in individs])


    def get_fitness(self):
        return selection.get_fitness(self)


    def hist_fitness(self):
        plt.hist(selection.get_fitness(self).values())



    def get_dom(self, chromosome, locus):
        return {locus: self.genomic_arch.dom[chromosome][locus]} 



    def get_env_var(self, chromosome, locus):
        return {locus: self.genomic_arch.env_var[chromosome][locus]} 





    def get_coords(self, individs = None):
        if individs <> None:
            return dict([(k, (ind.x, ind.y)) for k, ind in self.individs.items() if k in individs])
        else:
            return dict([(k, (ind.x, ind.y)) for k, ind in self.individs.items()])














    def show(self, land = None, scape_num = None, color = 'black', colorbar = True):
        if land <> None:
            if scape_num <> None  :
                land.scapes[scape_num].show(colorbar = colorbar)
            else:
                land.show(colorbar = colorbar)
        x = [(ind.x) for ind in self.individs.values()]
        y = [(ind.y) for ind in self.individs.values()]
        coords = np.array(self.get_coords().values()) - 0.5  #NOTE: subtract 0.5 to line up the points with the plt.imshow() grid of the land; imshow plots each pixel centered on its index, but the points then plot against those indices, so wind up shifted +-1.5
        mpl.pyplot.plot([n[0] for n in coords], [n[1] for n in coords], 'ko', scalex = False, scaley = False, color = color)




    def show_individs(self, individs, land = None, scape_num = None, color = 'black'):
        if land <> None and scape_num <> None:
            ax = land.scapes[scape_num].show()

        coords = dict([(k, (ind.x - 0.5, ind.y - 0.5)) for k, ind in self.individs.items() if k in individs]) #NOTE: subtract 0.5 to line up points with imshow grid; see note in the pop.show() definition for details
        for k, coord_pair in coords.items():
            ax = mpl.pyplot.plot(coord_pair[0], coord_pair[1], 'ko', scalex = False, scaley = False, color = color)
            #NOTE: perhaps worth figuring out how to label with the individual number!!







    def show_locus(self, chromosome, locus, land, scape_num = None):

        if scape_num <> None  :
            land.scapes[scape_num].show()
        else:
            land.show()

        genotypes = self.get_genotype(chromosome, locus)

        colors = ['#3C22B4', '#80A6FF', '#FFFFFF'] # COLORS TO MATCH LANDSCAPE PALETTE EXTREMES, BUT WITH HYBRID A MIX OF THE EXTREMES RATHER THAN THE YELLOW AT THE MIDDLE OF THE PALETTE, FOR NICER VIEWING: blue = [0,0], light blue = [0,1], white = [1,1]
        #colors = ['#ff4d4d', '#ac72ac', '#4d4dff'] # red = [0,0], purple = [0,1], blue = [1,1]

        for n, genotype in enumerate([0.0, 0.5, 1.0]):
            inds = [i for i, g in genotypes.items() if g[0] == genotype]
            coords = np.array([coord for coord in self.get_coords(inds).values()])

            mpl.pyplot.plot([coord[0] for coord in coords], [coord[1] for coord in coords], 'o', markersize = 11, scalex = False, scaley = False, color = colors[n])



    #method for plotting a population pyramid
    #NOTE: NEED TO FIX THIS SO THAT EACH HALF PLOTS ON OPPOSITE SIDES OF THE Y-AXIS
    def show_pyramid(self):
        plt.hist([ind.age for ind in self.individs.values() if ind.sex == 0], orientation = 'horizontal', color = 'pink', alpha = 0.6)
        plt.hist([ind.age for ind in self.individs.values() if ind.sex == 1], orientation = 'horizontal', color = 'skyblue', alpha = 0.6)








    def pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)











#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------



def create_population(genomic_arch, land, params):



    #grab necessary params from params dict

    N = params['N']

    dims = params['dims']

    size = params['size']

    T = params['T']






    assert dims.__class__.__name__ in ['tuple', 'list'], "dims should be expressed as a tuple or a list"
    individs = dict()
    for i in range(N):
        # use individual.create_individual to simulate individuals and add them to the population
        individs[i] = individual.create_individual(genomic_arch, dims)

    pop = Population(N = N, individs = individs, genomic_arch = genomic_arch, size = size, T = T)
    

    #get initial habitat values
    pop.query_habitat(land)

    return pop








#function for reading in a pickled pop
def load_pickled_pop(filename):
    import cPickle
    with open(filename, 'rb') as f:
        pop = cPickle.load(f)
    
    return pop
    


