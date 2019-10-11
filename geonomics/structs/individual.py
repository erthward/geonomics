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
from geonomics.structs.genome import _draw_genome
from geonomics.ops.movement import _do_dispersal
from geonomics.ops.selection import _calc_phenotype

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
        #individual's x-ploid genome
        self.g = new_genome
        #x and y coords
        self.x = float(x)
        self.y = float(y)
        #set individual's sex to that supplied, if supplied
        if sex:
            self.sex = sex
        #else, set it to a bernoulli r.v., p = 0.5;  0 --> female)
        else:
            self.sex = r.binomial(1, 0.5)
        self.age = age
        self.e = None
        self.z = []
        self.fit = None

        assert type(self.x) == float and self.x >= 0, ("invalid value "
                                "for x: %s, %s") % (str(self.x), type(self.x))
        assert type(self.y) == float and self.y >= 0, ("invalid value "
                                "for y: %s, %s") % (str(self.y), type(self.y))
        assert self.sex == None or self.sex in [0,1]
        assert type(self.age) == int


    #####################
    ### OTHER METHODS ###
    #####################

    #function to increment age by one
    def _set_age_stage(self):
        self.age += 1

    #sets the individual's position
    def _set_pos(self, x_pos, y_pos):
        self.x = x_pos
        self.y = y_pos

    #set the individual's current environmental values (attribute e)
    def _set_e(self, e):
        self.e = e

    #set the individual's phenotype (attribute z) for all traits
    def _set_z(self, genomic_architecture):
        self.z = [_calc_phenotype(self, genomic_architecture,
            trait) for trait in genomic_architecture.traits.values()]

    #set the individual's fitness
    def _set_fit(self, fit):
        self.fit = fit

    #set the individual's genome
    def _set_g(self, genome):
        self.g = genome


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

def _make_individual(idx, offspring=True, dim=None, genomic_architecture=None,
        new_genome = None, sex=None, parental_centerpoint = None,
        age=0, burn=False):

    """Create a new individual.

        If it is to have a genome, that can be created from either
        an instance of genome.Genomic_Architecture (i.e. for
        a newly simulated individual) or both the genomic architecture
        and a numpy.ndarray genome (e.g. for a new offspring).

        It will be assigned x- and y-coordinates, using the
        dim argument (i.e. the landscape dimensionality)
        if it is a newly simulated individual rather than offspring.

        It will be assigned a random sex, unless sex is provided.

        It will be assigned age 0, unless specifically fed otherwise.
        """
    #set the x,y location of the individual
    if offspring:
        #get a starting position from the parental_centerpoint
        x,y = _do_dispersal(parental_centerpoint)
    else:
        #randomly assign individual a valid starting location
        x,y = r.rand(2)*dim
        #clip to 0.01 under the dimensions, so that for landscapes even up to
        #~10000 on a side (this is bigger than I expect most users would
        #want to run) np.float32 can't return a number rounded up to
        #the dimension itself if an extremely high value is drawn
        x = np.clip(x, 0, dim[0]-0.001)
        y = np.clip(y, 0, dim[1]-0.001)

    #set the genome, if necessary
    if genomic_architecture is not None or new_genome is not None:
        #if not offspring (i.e. a new random individual), draw a new genome
        if not offspring:
            #if this is not for the burn-in, draw a proper genome
            if not burn:
                new_genome = _draw_genome(genomic_architecture)
            #otherwise, just a dummy genome
            else:
                new_genome = np.atleast_2d([0,0])

    #set the sex, if necessary
    if sex is None:
        #NOTE: For now sex randomly chosen at 50/50. Change if
        #decide to implement sex chroms, or spp.sex_ratio
        sex = r.binomial(1,0.5)

    return Individual(idx = idx, x = x, y = y, age = age,
                                        new_genome = new_genome, sex = sex)

