#!/usr/bin/python
#burn_in.py

import numpy as np
import matplotlib.pyplot as plt
from numpy import random as r
from scipy import interpolate
from scipy.spatial import cKDTree
from sklearn.preprocessing import normalize

import landscape
import mating








def burn_in(pop, land, params):

    #NOTE: Still not totally clear to me how to get this to work, but I want to use some iterative process to
    #get the model to empirically decide its own stable spatial population distribution based on its given
    #parameterization and start conditions, and to break the burn-in only once it has reached that state. I am
    #thinking that to do this I can use the previous timestep's population density as a local carrying capacity
    #for density-dependent births/deaths, and that by iterating the opposing forces of movement
    #to/through high-permeability areas and attrition in those areas because of density-dependence should reach
    #a dynamic equilibrium, which the model can determine through assessment of a metric or metrics




    #start the mean_dens array as just the starting pop density
    mean_dens = calc_pop_density(land, pop)


    #a while loop conditioned on burn == True, which will only switch to False once the stasis/stationarity criteria are determined to be met
    #while burn == True:
    for t in burn_T:

        #calculate pop density array, which will function as the carrying capacity (K) in each step
        K = calc_pop_density(land, pop)
        
        #calc local logistic growth rates

        
        #enact births and deaths, on basis of local logistic growth rates
            #NOTE: Probably should not be using selection here, as I want this to determine the initial
            #dynamic-equilibrium distribution on the basis of movement and density-dependence only


        #recalculate mean_dens (using the formula I worked out on the board for updating a mean)
        mean_dens = mean_dens + (1/(t+1))*(K-mean_dens)

        #assess stationarity of pop density
        
        #if determined stationary, based on chosen metric(s), 
            #set burn = False

        #else
            #repeat loop



