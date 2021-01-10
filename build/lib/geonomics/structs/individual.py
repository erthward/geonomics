#!/usr/bin/python
#individual.py

# flake8: noqa

'''
Defines the Individual class, with its associated methods and supporting
functions
'''

#geonomics imports
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
    """
    Representation of an individual of a given species.

    Multiple Individuals are collected as serial integer-keyed values within a
    Species dict. Those serial indices continually increment through model time
    as new Individuals (i.e. offspring) are produced, such that no two
    Individuals within the full history of a simulation will ever have the same
    index number.

    Attributes
    ----------

        NOTE: For more detail, see the documentation for the parameters
              that correspond to many of the following attributes.

        age:
            The Individual's age (expressed in timesteps since its birth
            timestep, when the individual was instantiated at age 0)

        e:
            The Individual's current environmental values, organized as a
            1d numpy array of length == len(Landscape), where the value stored
            at Individual.e[i] gives the Individual's current environmental
            value for Layer number i (i.e. Landscape[i])

        fit:
            The Individual's current fitness (calculated as a combination of
            the Individual's phenotypes for all traits, its current
            environmental values for all corresponding Landscape Layers, and
            any deleterious loci for which it has '1' alleles; see
            documentation for further details)

        g:
            The Individual's non-neutral genotypes,
            stored as an L_n x 2 numpy array, where L_n is the current number
            of non-neutral loci in the Individual's Species.
            Individual's only carry copies of their non-neutral genotypes
            (stored in this attribute). This is a computational optimization,
            as it allows fitness-based operations to be calculated
            quickly on the fly (using this attribute's data), while minimizing
            the memory required to stored the full (neutral and non-netural)
            genotypes of all current Individuals and their ancestors (i.e. the
            'spatial pedigree' stored in the succinct set of tskit tables).
            The successive rows in this array store the genotypes corresponding
            to the loci indicated by the successive locus numbers in
            Species.GenomicArchitecture.nonneut_loci.

        sex:
            The Individual's sex (0=female, 1=male; None, if the Species is
            unsexed)
        x:
            The Individual's current x (i.e. longitudinal) position
            (in continuous space, bounded within [0, Landscape.dims[0]])

        y:
            The Individual's current y (i.e. latitudinal) position
            (in continuous space, bounded within [0, Landscape.dims[1]])

        z:
            The Individual's phenotypes for each, organized as a
            1d numpy array of length == len(Species.traits), where the value
            stored at Individual.z[i] gives the Individual's fitness for
            Trait number i (i.e. Species.traits[i]).


    """
    def __init__(self, idx, x, y, age=0, new_genome=None, sex=None):
        self.idx = idx
        #individual's x-ploid genome (NOTE: make np.int8 to minimize mem)
        if new_genome is not None:
            self.g = np.int8(new_genome)
        else:
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

        # add attributes to hold the Individual's tskit Individuals id
        # and Nodes ids
        self._individuals_tab_id = None
        self._nodes_tab_ids = {}

        assert type(self.x) == float and self.x >= 0, ("invalid value "
                                "for x: %s, %s") % (str(self.x), type(self.x))
        assert type(self.y) == float and self.y >= 0, ("invalid value "
                                "for y: %s, %s") % (str(self.y), type(self.y))
        assert self.sex == None or self.sex in [0,1]
        assert type(self.age) == int


    #####################
    ### OTHER METHODS ###
    #####################

    # function to increment age by one
    def _set_age_stage(self):
        self.age += 1

    # set the individual's position
    def _set_pos(self, x_pos, y_pos):
        self.x = x_pos
        self.y = y_pos

    # set the individual's current environmental values (attribute e)
    def _set_e(self, e):
        self.e = e

    # set the individual's phenotype (attribute z) for all traits
    def _set_z(self, genomic_architecture):
        self.z = [_calc_phenotype(self, genomic_architecture,
            trait_num) for trait_num in genomic_architecture.traits]

    # set the individual's fitness
    def _set_fit(self, fit):
        self.fit = fit

    # set the individual's genome
    def _set_g(self, genome):
        self.g = genome

    # add a row of zeros to the individual's genome at the given index
    def _add_new_locus(self, idx):
        self.g = np.vstack((self.g[:idx, :],
                            [0, 0],
                            self.g[idx:, :]))

    # set the individual's tskit.TableCollection.nodes table's node ids
    # NOTE: node ids must be fed in 0-to-x order, where x is the species'
    # ploidy, because they will be assigned to the 0, 1, ..., x-keyed values
    # of the Individual._nodes_tab_ids dict according to that order
    def _set_nodes_tab_ids(self, *node_ids):
        self._nodes_tab_ids.update({i:id for i, id in enumerate(node_ids)})

    # get the individual's x,y location as a tuple
    def _get_loc(self):
        return(self.x, self.y)


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

def _make_individual(idx, offspring=True, dim=None, genomic_architecture=None,
        new_genome=None, sex=None, parental_centerpoint=None,
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

    #set the sex, if necessary
    if sex is None:
        #NOTE: For now sex randomly chosen at 50/50. Change if
        #decide to implement sex chroms, or spp.sex_ratio
        sex = r.binomial(1,0.5)

    individual = Individual(idx=idx, x=x, y=y, age=age,
                            new_genome=new_genome, sex=sex)
    return individual

