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

            


land = landscape.make_land(params = params)
com = community.make_community(land = land, params = params, burn = True)
pop = com[0]

