#!/usr/bin/python
#burn_in.py

'''Functions for testing burn-in stationarity.'''

import numpy as np
import matplotlib.pyplot as plt
from numpy import random as r
from scipy import interpolate
from scipy.spatial import cKDTree
from scipy.stats import ttest_rel as tt
from statsmodels.tsa.stattools import adfuller as adf
from sklearn.preprocessing import normalize

from structs import landscape
from ops import mating


#TODO: Add more and/or different tests? This is a basic test so far, based
#only on gross population size (not spatial distribution)


def test_adf_threshold(pop, num_timesteps_back, p_val=0.05):
    return(adf(pop.Nt[-num_timesteps_back:])[1] < p_val)


def test_t_threshold(pop, num_timesteps_back, p_val=0.05):
    num_timesteps_back += num_timesteps_back % 2
    return(tt(pop.Nt[int(-num_timesteps_back): int(-num_timesteps_back/2)], pop.Nt[int(-num_timesteps_back/2):])[1] >0.05)




