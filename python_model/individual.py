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
    def __init__(self, G, x, y, sex = None, age=0):

        self.genome = G             #individual's x-ploid genome

        self.x = float(x)           #x coord
        self.y = float(y)           #y coord

        if sex:
            self.sex = sex          #set individual's sex to that supplied, if supplied
        else:
            self.sex = r.binomial(1, 0.5) #else, set it to a bernoulli r.v., p = 0.5;  0 --> female)

        self.age = age              #age

        self.habitat = None



        assert self.genome.__class__.__name__ == 'Genome', "An individual's genome must be an instance of the genome.Genome class."
        assert type(self.x) == float and self.x >= 0, "invalid value for x: %s, %s" % (str(self.x), type(self.x))
        assert type(self.y) == float and self.y >= 0, "invalid value for y: %s, %s" % (str(self.y), type(self.y))
        assert self.sex == None or self.sex in [0,1]
        assert type(self.age) == int





    #####################
    ### OTHER METHODS ###
    #####################


    #set movement.move as a method
    def move(self, land, params):
        self.x, self.y = movement.move(self, land, params)
        

    
    #function to increment age by one
    def increment_age_stage(self):
        self.age = self.age + 1


    #function to extract the habitat value for an individual's current raster cell
    def query_habitat(self, land):
        if self.habitat == None:
            self.habitat = [np.nan] * land.num_rasters
        for n in land.scapes.keys():
            j, i = [int(coord) for coord in np.floor([self.x, self.y])]  #NOTE: j,i NOT i,j becuase the numpy array will be indexed i then j, i.e. vertical, then horizontal, which corresponds to y then x!!  
            self.habitat[n] = land.scapes[n].raster[i,j]






#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------


def create_individual(genomic_arch, dims=None, genomic_content = None, ploidy = None, parental_centerpoint = None, sex = None, age=0):
    '''Create a new individual from:
            - either an instance of genome.Genomic_Architecture (i.e. for newly simulated individual) or both the genomic architecture and an instance of genome.Genome (e.g. for new offspring) (one of the two must be provided),
            - x and y-coordinates
            - sex
            - age.
            '''

    #LOOP FOR SIMULATION OF NEW INDIVIDUALS FOR STARTING POPULATION
    if genomic_content == None:
        assert dims <> None, "landscape dims required to simulate a new individual"
        #use genome.sim_genome and genomic_arch variable to simulate individual's genome
        G = genome.sim_genome(genomic_arch)
         
        #randomly assign individual a valid starting location
        x,y = r.rand(2)*dims
        return Individual(G, x, y, sex = sex, age = age)


    #LOOP FOR CREATION OF NEW OFFSPRING INDIVIDUALS
    else:
        assert parental_centerpoint <> None, "parental_centerpoint needed to create new offspring"
        assert parental_centerpoint.__class__.__name__ in ['tuple', 'list'], "parental_centerpoint should be a tuple or a list"
        assert parental_centerpoint[0] >= 0 and parental_centerpoint[1] >= 0, "parental_centerpoint coordinates must be within landscape, but %s was provided" % str(parental_centerpoint)
        assert ploidy <> None, "ploidy needed to create new genome from genomic content, for assertions in genome.Genome class __init__ function"


        x,y = dispersal.disperse(parental_centerpoint) #NOTE: needs to be written!


        sex = r.binomail(1,0.5)  #NOTE: For now, sex randomly chosen at 50/50. Change if later decide to implement sex chroms!!!


        return Individual(G, x, y, sex = sex, age = 0)


