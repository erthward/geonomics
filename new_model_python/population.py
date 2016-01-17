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


class Population:
    def __init__(self, N, individs, birth_rate, death_rate, genomic_arch, T):
        self.N = [N]                                        #list to record population size at each step (with starting size as first entry)
        self.individs = individs                            #dict of all individuals, as instances of individual.Individual class
        if type(birth_rate) == float:
            self.birth_rate = [birth_rate] * T              #birth rate, expressed as vector of T values, where T is total model runtime; can be input as constant value, or as vector varying over model
        elif type(birth_rate) == list and len(birth_rate) == T:
            self.birth_rate = birth_rate
        else:
            assert type(birth_rate) == list, "birth_rate should be expressed either as a floating point value between 0 and 1, or as a list of such values"
            assert len(birth_rate) == T, "if expressing birth_rate as a list, list length must equal total model runtime"
        if type(death_rate) == float:                       #death rate, handled same as birth rate
            self.death_rate = [death_rate] * T              
        elif type(death_rate) == list and len(death_rate) == T:
            self.death_rate = death_rate
        else:
            assert type(death_rate) == list, "death_rate should be expressed either as a floating point value between 0 and 1, or as a list of such values"
            assert len(death_rate) == T, "if expressing death_rate as a list, list length must equal total model runtime"

        self.death_rate = death_rate   
        self.genomic_arch = genomic_arch
        #Add other demographic parameters?? Other stuff in general??


        assert type(N) == int, "N must be an integer"
        assert type(self.individs) == dict, "self.individs must be a dictionary"
        assert list(set([i.__class__.__name__ for i in self.individs.values()])) == ['Individual'], "self.individs must be a dictionary containing only instances of the individual.Individual class"
        assert self.genomic_arch.__class__.__name__ == 'Genomic_Architecture', "self.genomic_arch must be an instance of the genome.Genomic_Architecture class"

    def show(self, land = None, color = 'black'):
        if land:
            land.show()
        x = [(ind.x) for ind in self.individs.values()]
        y = [(ind.y) for ind in self.individs.values()]
        mpl.pyplot.plot(x,y, 'ko', scalex = False, scaley = False, color = color)


    #method to move all individuals simultaneously
    def move(self, land):
        [ind.move(land) for ind in self.individs.values()];


    #set mating.find_mates() as method
    def find_mates(self, mating_radius, sex = True, repro_age = None):
        return mating.find_mates(self, mating_radius, sex = sex, repro_age = repro_age)


    #function for executing mating for a population
    def mate(self, land, mating_radius, mu_dispersal, sigma_dispersal, sex = True, repro_age = None):

        mating_pairs = self.find_mates(mating_radius, sex = sex, repro_age = repro_age)

        print '\n\nMATING PAIRS:\n\n'
        print mating_pairs

        for pair in mating_pairs:

            zygote = mating.mate(self, pair, self.genomic_arch)

            parent_centroid_x = np.mean((self.individs[pair[0]].x, self.individs[pair[1]].x))
            parent_centroid_y = np.mean((self.individs[pair[0]].y, self.individs[pair[1]].y))
            offspring_x, offspring_y = dispersal.disperse(land, parent_centroid_x, parent_centroid_y, mu_dispersal, sigma_dispersal)

            sex = r.binomial(1, 0.5)

            age = 0

            offspring_num = max(self.individs.keys()) + 1

            self.individs[offspring_num] = individual.Individual(zygote, offspring_x, offspring_y, sex, age)


            



    #method to increment all population's age by one
    def birthday(self):
        [ind.birthday() for ind in self.individs.values()];



    def pickle(self, filename):
        import cPickle as pickle
        with open(filename, 'wb') as f:
            pickle.dump(self, f)




def create_population(N, genomic_arch, dims, birth_rate, death_rate, T, sex = False):
    assert type(dims) in [tuple, list], "dims should be expressed as a tuple or a list"
    individs = dict()
    for i in range(N):
        # use individual.create_individual to simulate individuals and add them to the population
        individs[i] = individual.create_individual(genomic_arch, dims)

    return Population(N, individs, birth_rate, death_rate, genomic_arch, T)

