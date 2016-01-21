#!/usr/bin/python 
#population.jl

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

import numpy as np
from numpy import random as r
import matplotlib as mpl

#------------------------------------
# CLASSES ---------------------------
#------------------------------------


class Population:



    '''Population(self, N, individs, genomic_arch, size, birth_rate = None, death_rate = None, T = None)'''


     
    def __init__(self, N, individs, genomic_arch, size, birth_rate = None, death_rate = None, T = None): #NOTE: PROBABLY DELETE THE LAST THREE ARGS!!
        self.N = [N]                                        #list to record population size at each step (with starting size as first entry)
        self.individs = individs                            #dict of all individuals, as instances of individual.Individual class





##########################################################
#NOTE: BEGIN SECTION TO PROBABLY DELETE!!

        if birth_rate <> None:
            if type(birth_rate) == float:
                self.birth_rate = [birth_rate] * T              #birth rate, expressed as vector of T values, where T is total model runtime; can be input as constant value, or as vector varying over model
            elif type(birth_rate) == list and len(birth_rate) == T:
                self.birth_rate = birth_rate
            else:
                assert type(birth_rate) == list, "birth_rate should be expressed either as a floating point value between 0 and 1, or as a list of such values"
                assert len(birth_rate) == T, "if expressing birth_rate as a list, list length must equal total model runtime"
        if death_rate <> None:
            if type(death_rate) == float:                       #death rate, handled same as birth rate
                self.death_rate = [death_rate] * T              
            elif type(death_rate) == list and len(death_rate) == T:
                self.death_rate = death_rate
            else:
                assert type(death_rate) == list, "death_rate should be expressed either as a floating point value between 0 and 1, or as a list of such values"
                assert len(death_rate) == T, "if expressing death_rate as a list, list length must equal total model runtime"


#END SECTION TO PROBABLY DELETE
##########################################################





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



    #method to extract habitat values at each individual's coordinates for each land.scapes raster
    def query_habitat(self, land):
        [ind.query_habitat(land) for ind in self.individs.values()];
 

    #method to move all individuals simultaneously
    def move(self, land, mu_direction = 0, kappa_direction = 0, mu_distance = 4, sigma_distance = 0.1, resist_surf = None):
        [ind.move(land, mu_direction, kappa_direction, mu_distance, sigma_distance, resist_surf) for ind in self.individs.values()];
        self.query_habitat(land)




    #set mating.find_mates() as method
    def court(self, mating_radius, sex = True, repro_age = None):
        return mating.find_mates(self, mating_radius, sex = sex, repro_age = repro_age)




    #function for executing mating for a population
    def mate(self, land, mating_radius, mu_dispersal, sigma_dispersal, sex = True, repro_age = None):

        mating_pairs = self.court(mating_radius, sex = sex, repro_age = repro_age)

        for pair in mating_pairs:

            parent_centroid_x = np.mean((self.individs[pair[0]].x, self.individs[pair[1]].x))
            parent_centroid_y = np.mean((self.individs[pair[0]].y, self.individs[pair[1]].y))

            zygotes = mating.mate(self, pair, self.genomic_arch)

            
            for zygote in zygotes:
                
                offspring_x, offspring_y = dispersal.disperse(land, parent_centroid_x, parent_centroid_y, mu_dispersal, sigma_dispersal)


                sex = r.binomial(1, 0.5)

                age = 0

                offspring_key = max(self.individs.keys()) + 1

                self.individs[offspring_key] = individual.Individual(zygote, offspring_x, offspring_y, sex, age)


        #sample all individuals' habitat values, to initiate for offspring
        self.query_habitat(land)




                


            



    #method to increment all population's age by one
    def birthday(self):
        [ind.birthday() for ind in self.individs.values()];


    

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
 

   


   
 


    def show(self, land = None, scape_num = None, color = 'black'):
        if land <> None:
            if scape_num <> None  :
                land.scapes[scape_num].show()
            else:
                land.show()
        x = [(ind.x) for ind in self.individs.values()]
        y = [(ind.y) for ind in self.individs.values()]
        coords = np.array(self.get_coords().values()) - 0.5  #NOTE: subtract 0.5 to line up the points with the plt.imshow() grid of the land; imshow plots each pixel centered on its index, but the points then plot against those indices, so wind up shifted +-1.5
        print coords
        mpl.pyplot.plot([n[0] for n in coords], [n[1] for n in coords], 'ko', scalex = False, scaley = False, color = color)




    def show_individs(self, individs, land = None, scape_num = None, color = 'black'):
        if land <> None and scape_num <> None:
            ax = land.scapes[scape_num].show()

        coords = dict([(k, (ind.x - 0.5, ind.y - 0.5)) for k, ind in self.individs.items() if k in individs]) #NOTE: subtract 0.5 to line up points with imshow grid; see note in the pop.show() definition for details
        for k, coord_pair in coords.items():
            ax = mpl.pyplot.plot(coord_pair[0], coord_pair[1], 'ko', scalex = False, scaley = False, color = color)
            #NOTE: perhaps worth figuring out how to label with the individual number!!




    def get_habitat(self, individs = None):
        if individs <> None:
            return dict([(k, ind.habitat) for k, ind in self.individs.items() if k in individs])
        else:
            return dict([(k, ind.habitat) for k, ind in self.individs.items()])


    

    def get_genotype(self, locus, chromosome_num, format = 'mean', individs = None):

        dom_funcs = { 0 : np.mean,                          #codominance
                      1 : lambda x: np.ceil(np.mean(x)),    #dominant (with respect to relationship between allele 1 and habitat values -> 1
                      2 : lambda x: np.floor(np.mean(x))    #recessive (with respect to relationship between allele 1 and habitat values -> 1
                    }


        if individs == None:
            individs = range(len(self.genomic_arch.s[chromosome_num]))

        if format == 'biallelic':
            return dict([(i, self.individs[i].genome.genome[chromosome_num][locus, :]) for i in self.individs.keys() if i in individs]) 

        elif format == 'mean':
            return dict([(i, [dom_funcs[self.genomic_arch.dom[chromosome_num][locus]](self.individs[i].genome.genome[chromosome_num][locus,:]) ]) for i in self.individs.keys() if i in individs])




    def get_dom(self, chromosome_num, locus):
        return {locus: self.genomic_arch.dom[chromosome_num][locus]} 






    def get_coords(self, individs = None):
        if individs <> None:
            return dict([(k, (ind.x, ind.y)) for k, ind in self.individs.items() if k in individs])
        else:
            return dict([(k, (ind.x, ind.y)) for k, ind in self.individs.items()])






    def pickle(self, filename):
        import cPickle
        with open(filename, 'wb') as f:
            cPickle.dump(self, f)



#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------



def create_population(N, genomic_arch, dims, size, land, birth_rate = None, death_rate = None, T = None):
    assert dims.__class__.__name__ in ['tuple', 'list'], "dims should be expressed as a tuple or a list"
    individs = dict()
    for i in range(N):
        # use individual.create_individual to simulate individuals and add them to the population
        individs[i] = individual.create_individual(genomic_arch, dims)

    pop = Population(N = N, individs = individs, genomic_arch = genomic_arch, size = size, T = T)
    

    #get initial habitat values
    pop.query_habitat(land)

    return pop



def load_pickled_pop(filename):
    import cPickle
    with open(filename, 'rb') as f:
        pop = cPickle.load(f)
    
    return pop
    


