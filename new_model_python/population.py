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
        assert self.genomic_arch.__class__.__name__ == 'Genomic_architecture', "self.genomic_arch must be an instance of the genome.Genomic_Architecture class"

    def show(self):
        x = [(ind.x) for ind in self.individs.values()]
        y = [(ind.y) for ind in self.individs.values()]
        mpl.pyplot.plot(x,y, 'ko', scalex = False, scaley = False)




def create_population(N, genomic_arch, dims, birth_rate, death_rate, T, sex = False):
    assert type(dims) in [tuple, list], "dims should be expressed as a tuple or a list"
    individs = dict()
    for i in range(N):
        # use individual.create_individual to simulate individuals and add them to the population
        individs[i] = individual.create_individual(genomic_arch, dims)

    return Population(N, individs, birth_rate, death_rate, genomic_arch, T)

