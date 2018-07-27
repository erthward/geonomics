#!/usr/bin/python

#imports_and_reloads.py

'''Import and reload all modules, to call when working on the command line.'''

from importlib import reload

from utils import viz, io, spatial as spt
from structs import landscape, genome, individual, population, community
from ops import movement, demography, mating, selection, mutation, change
from sims import burn_in, stats, data
from support import param_exploration as pe

import numpy as np
from numpy import random as r
import random
import matplotlib as mpl
import os

for pkg in [viz, io, spt, landscape, genome, individual, population, community, movement, demography, mating, selection, mutation, burn_in, stats, data, change, pe]:
    reload(pkg)
