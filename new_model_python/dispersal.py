#!/usr/bin/python
#dispersal.py


'''
##########################################

Module name:                dispersal.py

Module contains:
                            - function for simulating the dispersal of offspring, according to input parameters
                            - associated functions


Author:                     Drew Ellison Hart
Email:                      drew.hart@berkeley.edu
Github:                     URL
Start date:                 12-28-15
Documentation:              URL


##########################################
'''

import numpy as np
from numpy.random import vonmises, lognormal



#------------------------------------
# CLASSES ---------------------------
#------------------------------------


#--------------------------------------
# FUNCTIONS ---------------------------
#--------------------------------------



def disperse(land, parent_centroid_x, parent_centroid_y, mu_dispersal, sigma_dispersal, mu_dir = 0, kappa_dir = 0): 


    within_landscape = False
    while within_landscape == False:

        direction = vonmises(mu_dir, kappa_dir)   #NOTE: For now, dispersal random and equally probable in all directions 
        distance = lognormal(mu_dispersal, sigma_dispersal)

        offspring_x = parent_centroid_x + np.cos(direction)*distance
        offspring_y = parent_centroid_y + np.sin(direction)*distance
        within_landscape = (offspring_x > 0 and offspring_x < land.dims[0]) and (offspring_y > 0 and offspring_y < land.dims[1])

    return (offspring_x, offspring_y)


   




