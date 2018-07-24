#/usr/bin/python
#PRACTICE_create_land_genomic_arch_pop.py

'''Creates a landscape object, a genomic architecture object, and population object, using the current params
dictionary.'''


import numpy.random as r
from structs import landscape, genome, population


#set seed, if requested
if params['seed']['set_seed']: 
    random.seed(params['seed']['seed_num'])
    r.seed(params['seed']['seed_num'])

            


land = landscape.make_scape_stack(params = params)
arch = genome.make_genomic_arch(params = params, land = land)
com = population.make_community(genomic_arch = arch, land = land, params = params)

pop = com[0]

