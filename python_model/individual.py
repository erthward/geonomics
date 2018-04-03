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


import genome
import movement


import numpy as np
import numpy.random as r


#------------------------------------
# CLASSES ---------------------------
#------------------------------------


class Individual:
    def __init__(self, new_genome, x, y, sex = None, age=0):

        self.genome = new_genome             #individual's x-ploid genome

        self.x = float(x)           #x coord
        self.y = float(y)           #y coord

        if sex:
            self.sex = sex          #set individual's sex to that supplied, if supplied
        else:
            self.sex = r.binomial(1, 0.5) #else, set it to a bernoulli r.v., p = 0.5;  0 --> female)

        self.age = age              #age

        self.habitat = None



        assert type(self.genome) == np.ndarray, "An individual's genome must be an instance of numpy.ndarray."
        assert type(self.x) == float and self.x >= 0, "invalid value for x: %s, %s" % (str(self.x), type(self.x))
        assert type(self.y) == float and self.y >= 0, "invalid value for y: %s, %s" % (str(self.y), type(self.y))
        assert self.sex == None or self.sex in [0,1]
        assert type(self.age) == int





    #####################
    ### OTHER METHODS ###
    #####################


    #set movement.move as a method
    def move(self, land, params):
        movement.move(self, land, params)


    #function to increment age by one
    def increment_age_stage(self):
        self.age += 1

    def set_pos(self, pos_x, pos_y):
        self.x = pos_x
        self.y = pos_y




#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------


def create_individual(genomic_arch, dims=None, new_genome = None, ploidy = None, parental_centerpoint = None, sex = None, age=0, burn =False):
    """Create a new individual from:
            - either an instance of genome.Genomic_Architecture (i.e. for newly simulated individual) or both
              the genomic architecture and a numpy.ndarray genome (e.g. for new offspring) (one of the two must be provided),
            - x and y-coordinates
            - sex
            - age.
            """

    #LOOP FOR SIMULATION OF NEW INDIVIDUALS FOR STARTING POPULATION

    if new_genome != None:
        assert parental_centerpoint != None, "parental_centerpoint needed to create new offspring"
        assert parental_centerpoint.__class__.__name__ in ['tuple', 'list'], "parental_centerpoint should be a tuple or a list"
        assert parental_centerpoint[0] >= 0 and parental_centerpoint[1] >= 0, "parental_centerpoint coordinates must be within landscape, but %s was provided" % str(parental_centerpoint)
        assert ploidy != None, "ploidy needed to create new genome from genomic content"


        x,y = dispersal.disperse(parental_centerpoint) #NOTE: needs to be written!


        sex = r.binomial(1,0.5)  #NOTE: For now, sex randomly chosen at 50/50. Change if later decide to implement sex chroms!!!


        return Individual(new_genome, x, y, sex = sex, age = 0)


    elif new_genome == None:
        assert dims != None, "landscape dims required to simulate a new individual without reproduction"

        #randomly assign individual a valid starting location
        x,y = r.rand(2)*dims

        if burn == False:

            #use genome.sim_genome and genomic_arch variable to simulate individual's genome
            new_genome = genome.sim_genome(genomic_arch)

            return Individual(new_genome, x, y, sex = sex, age = age)

        elif burn == True:

            return(Individual(new_genome = np.array([0]), x=x, y=y, sex = sex, age = age))


