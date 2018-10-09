#!/usr/bin/python
#individual.py

'''
##########################################

Module name:              individual


Module contains:
                          - definition of the Individual type
                          - function for creation of a simulated individual
                          - associated functions


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               12-28-15
Documentation:            URL


##########################################
'''

#geonomics imports
from structs import genome
from ops import movement, selection

#other imports
import numpy as np
import numpy.random as r


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

class Individual:
    def __init__(self, idx, x, y, age=0, new_genome=None, sex=None):

        self.idx = idx

        self.genome = new_genome             #individual's x-ploid genome

        self.x = float(x)           #x coord
        self.y = float(y)           #y coord

        if sex:
            self.sex = sex          #set individual's sex to that supplied, if supplied
        else:
            self.sex = r.binomial(1, 0.5) #else, set it to a bernoulli r.v., p = 0.5;  0 --> female)

        self.age = age              #age

        self.habitat = None

        self.phenotype = []

        self.fitness = None



        #assert type(self.genome) == np.ndarray or self.genome is None, "An individual's genome must either be None or an instance of numpy.ndarray."
        assert type(self.x) == float and self.x >= 0, "invalid value for x: %s, %s" % (str(self.x), type(self.x))
        assert type(self.y) == float and self.y >= 0, "invalid value for y: %s, %s" % (str(self.y), type(self.y))
        assert self.sex == None or self.sex in [0,1]
        assert type(self.age) == int





    #####################
    ### OTHER METHODS ###
    #####################


    #function to increment age by one
    def set_age_stage(self):
        self.age += 1

    #sets the individual's position
    def set_pos(self, x_pos, y_pos):
        self.x = x_pos
        self.y = y_pos

    #set the individual's habitat
    def set_habitat(self, hab):
        self.habitat = hab

    #set the individual's phenotype for all traits
    def set_phenotype(self, genomic_architecture):
        self.phenotype = [selection.calc_phenotype(self, genomic_architecture, 
            trait) for trait in genomic_architecture.traits.values()]

    #set the individual's fitness
    def set_fitness(self, fit):
        self.fitness = fit

    #set the individual's genome
    def set_genome(self, genome):
        self.genome = genome


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

def make_individual(idx, offspring=True, dim=None, genomic_architecture=None, new_genome = None, sex=None, parental_centerpoint = None, age=0, burn=False):

    """Create a new individual.

        If it is to have a genome, that can be created from either an instance of 
        genome.Genomic_Architecture (i.e. for a newly simulated individual) or 
        both the genomic architecture and a numpy.ndarray genome (e.g. for a new 
        offspring).
    
        It will be assigned x- and y-coordinates, using the dim argument (i.e. the
        landscape dimensionality) if it is a newly simulated individual rather 
        than offspring.

        It will be assigned a random sex, unless sex is provided.

        It will be assigned age 0, unless specifically fed otherwise.
        """
    #set the x,y location of the individual
    if offspring:
        #TODO: probably remove these assert statements; don't see any reason I need to keep running them, so it's just unnecessary comuptation during the model
        #assert parental_centerpoint != None, "parental_centerpoint needed to create new offspring"
        #assert parental_centerpoint.__class__.__name__ in ['tuple', 'list'], "parental_centerpoint should be a tuple or a list"
        #assert parental_centerpoint[0] >= 0 and parental_centerpoint[1] >= 0, "parental_centerpoint coordinates must be within landscape, but %s was provided" % str(parental_centerpoint)
        #get a starting position from the parental_centerpoint
        x,y = dispersal.disperse(parental_centerpoint) 
    else:
        #randomly assign individual a valid starting location
        x,y = r.rand(2)*dim

    #set the genome, if necessary
    if genomic_architecture is not None or new_genome is not None:
        #if not offspring (i.e. a new random individual), draw a new genome
        if not offspring:
            #if this is not for the burn-in, draw a proper genome
            if not burn:
                new_genome = genome.draw_genome(genomic_architecture)
            #otherwise, just a dummy genome
            else:
                new_genome = np.atleast_2d([0,0])

    #set the sex, if necessary
    if sex is None:
        #NOTE: For now sex randomly chosen at 50/50. Change if decide to implement sex chroms, or pop.sex_ratio
        sex = r.binomial(1,0.5)  

    return Individual(idx = idx, x = x, y = y, age = age, new_genome = new_genome, sex = sex)

