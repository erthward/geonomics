#!/usr/bin/python
#burn_in.py

'''Helper functions for running and breaking burn-in.'''

import numpy as np
import matplotlib.pyplot as plt
from numpy import random as r
from scipy import interpolate
from scipy.spatial import cKDTree
from scipy.stats import ttest_rel as tt
from statsmodels.tsa.stattools import adfuller as adf
from sklearn.preprocessing import normalize

import landscape
import mating



def adf_threshold_test(pop, num_timesteps_back, p_val):
    return(adf(pop.Nt[-num_timesteps_back:])[1] < p_val)


def tt_threshold_test(pop, num_timesteps_back, p_val):
    num_timesteps_back += num_timesteps_back % 2
    return(tt(pop.Nt[int(-num_timesteps_back): int(-num_timesteps_back/2)], pop.Nt[int(-num_timesteps_back/2):])[1] >0.05)




