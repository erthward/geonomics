#!/usr/bin/python
#demography.py

'''
##########################################

Module name:          ops.demography


Module contains:
                      - Functions to control demography/population dynamics
                      - a DebugPlotter class, which will create a series
                        of plots to help debug the pop_dynamics function
                        (and general model behavior) if that function is
                        called with debug = True
 

Author:               Drew Ellison Hart
Email:                drew.hart@berkeley.edu
Github:               URL
Start date:           04-28-17
Documentation:        URL


##########################################
'''

#geonomics imports
from structs import landscape
from ops import mating, selection

#other imports
import numpy as np
import matplotlib.pyplot as plt
from numpy import random as r
from scipy import interpolate
from sklearn.preprocessing import normalize
from collections import OrderedDict as OD
from operator import itemgetter


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

class _DebugPlotter:
    def __init__(self, land, pop, timestep):
        self.subplot = 241
        self.timestep = timestep
        self.pop = pop
        self.land = land

    def _next_plot(self, varname, var):
        if self.subplot <= 248:
            if not plt.get_fignums():
                fig = plt.figure()
                plt.suptitle('TIMESTEP: %i' % self.timestep)
            else:
                fig = plt.gcf()
            plt.subplot(self.subplot)
            plt.title(varname)
            self.pop.plot(self.land, 1, size = 8, hide_land = True)
            plt.imshow(var, interpolation = 'nearest', cmap = 'terrain')
            plt.colorbar()
            self.subplot += 1
        else:
            pass


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

def _calc_n_pairs(pairs, pop):
    #if there are pairs
    if pairs.size > 0:
        #get their coordinates
        pairs_coords = pop._get_coords(individs = pairs.flatten())
        #take just the x and y coordinates from pairs_coords, reshape them to
        #match pairs.shape (i.e.nx2, where n = len(pop)), and then mean across
        #axis 1 (to get the pair's centroid
        #x and y coordinates)
        p_x = pairs_coords[:,0].reshape(pairs.shape).mean(axis = 1)
        p_y = pairs_coords[:,1].reshape(pairs.shape).mean(axis = 1)

        #use the pop._dens_grids to calculate a raster of the expected number
        #of pairs
        #NOTE: because the window_width of that object is much wider than 
        #dispersal_mu in the parameterizations I've typically been running, 
        #this produces much wider areas of high values in the N_b, N_d, and
        #hence d rasters than the old method did (which had its window_width
        #set to dispersal_mu, leading of course to huge lags in the old 
        #calc_density algorithm) need to think about whether I need to attend
        #to this, and if so, how
        n_pairs = np.clip(pop._dens_grids._calc_density(p_x, p_y), a_min = 0,
                                                                a_max = None)

        #remove NaNs
        n_pairs[np.isnan(n_pairs)] = 0

    #otherwise just return a 0 raster
    else:
        n_pairs = np.zeros(pop._land_dim)

    return n_pairs


#the logistic eqxn
def _calc_logistic_growth(R, N, K):
    dNdt = R*(1-(N/K))*N
    return dNdt


def _calc_logistic_soln(x0, R, t):
        return 1/(1+((1/x0)-1)*np.e**(-1*R*t))


def _calc_dNdt(land, R, N, K, pop_growth_eq = 'logistic'):
    with np.errstate(divide='ignore', invalid='ignore'):
        dim = land.dim
        #if using the logistic equation for pop growth
        #NOTE: For now this is the only option and the default, but could 
        #easily offer other equations later if desired
        if pop_growth_eq == 'logistic':
            #use logistic eqxn, with pop intrinsic growth rate (pop.R)
            #to generate current growth-rate raster
            dNdt = _calc_logistic_growth(R, N, K)
    #NOTE: The next line used to replace each cell in dNdt where K<1
    #with that cell's -N value. But it was really slow to run, and at those
    #cells where K ~= 0 it really won't make a difference how negative the
    #dNdt values there are, so the following line just makes this much simpler.
    dNdt = np.clip(dNdt, a_min = -1*N.max(), a_max = None)
    #NOTE: DEH 06-07-18: I haven't thought it through well enough to guarantee
    #that this is the best way to handle NaNs and infs, but they at least result
    #from division by 0 zero cells in the pop.K raster, so for now
    #correcting them to the minimum val
    dNdt[np.isnan(dNdt)] = -1*N.max()
    dNdt[np.isinf(dNdt)] = -1*N.max()
    return dNdt


def _calc_N_b(b, n_births_distr_lambda, n_pairs):
    #Use n_pairs, pop.b (i.e. probability that a pair gives birth),
    #and pop.n_births_distr_lambda (i.e. expected number of births
    #per pair) to calculate the N_b raster (expected number of births in each
    #cell) according to N_b = b*E[N_b/pair]*n_pairs where E[N_b/pair] = 
    #n_births_distr_lambda (b/c births Poission distributed)
    #NOTE: Would be good to in the future to use neg_binom dist to allow
    #births to be modeled as Poisson (using correct parameterization of
    #neg_binom to approximate that) or overdispersed (to be parmeterized
    #with field data)

    #NOTE: Does it matter that I'm not attempting to account for the fact
    #that pairs will disperse out of the cell? Almost certainly a larger effect
    #with increasing dispersal-kernel distance. But nonetheless, would it even
    #make sense (and even BE TRACTABLE) to calculate it any other way? My
    #intuition is NO, but worth considering further later on...

    #NOTE: Currently, birth rate is intrinsic to the species and cannot be
    #made a function of density or spatial fitness. Might be worth adding
    #this at some point in the future though.
    N_b = b * n_births_distr_lambda * n_pairs
    return N_b


def _calc_Nd(N_b, dNdt):
    #Use N_b - dN/dt to calculate N_d raster (expected number of deaths in
    #each cell)
    N_d = N_b - dNdt
    return N_d


def _calc_d(N_d, N, d_min, d_max):
    #Use N_d/N to calculate d raster (probability of death of individuals
    #in each cell) Calculate d as simply N_d/N, then control for some
    #numerical artefacts that an inappropriate
    d = N_d/N
    #fix infinties and NaNs and negatives if they arise
    #(NOTE: infinity occurs where N == 0)
    d[np.isinf(d)] = 0
    d[np.isnan(d)] = 0
    d[d<0] = 0
    #constrain to the min and max d values
    d[d<d_min] = d_min
    d[d>d_max] = d_max
    return d


def _do_mortality(land, pop, death_probs):
    deaths = np.array([*pop])[np.bool8(r.binomial(n = 1, p = death_probs))]
    if len(deaths) > 0:
        ig = itemgetter(*deaths)
        [pop.pop(ind) for ind in deaths];
    return len(deaths)


def _do_pop_dynamics(land, pop, with_selection = True, burn = False,
                births_before_deaths = False, asserts = False, debug = False):
    '''Generalized function for implementation population dynamics.
    Will carry out one round of mating and death, according to parameterization
    laid out in params dict (which were grabbed as Population attributes).

       If with_selection == False, only density-dependent death will occur.

       If burn == True, selection will be coerced to False, the Population.K
       attribute will be updated, stationarity metrics will be assessed, and
       the function will return the decision it has reached regarding
       stationarity.
    '''

    #NOTE: Would be good to write in option to ONLY mate or ONLY die when the
    #function is run (so that a pop could be run with mating every x timesteps
    #but mortality each timestep, for example)

    #create the debug-plotter, id debug == True
    if debug:
        timestep = -99999
        dp = _DebugPlotter(land, pop, timestep)

    #find mating pairs
    pairs = pop._find_mating_pairs()

    #calc num_pairs raster (use the calc_pop_density function on the
    #centroids of the mating pairs)
    n_pairs = _calc_n_pairs(pairs = pairs, pop = pop)
    #run checks on n_pairs
    if asserts:
        assert n_pairs.min() >= 0, 'n_pairs.min() == %0.2f' %(n_pairs.min())
        assert True not in np.isnan(n_pairs)
        assert True not in np.isinf(n_pairs)
    #add debug plot
    if debug:
        dp._next_plot('n_pairs', n_pairs)

    #if births should happen before (and thus be included in the calculation of)
    #deaths, then mate and disperse babies now
    if births_before_deaths:
        #Feed the land and mating pairs to pop.do_mating, to produce and
        #disperse zygotes
        pop._do_mating(pairs, burn)

    #calc N raster, set it as pop.N, then get it  
    pop._calc_density(set_N = True)
    N = pop.N
    #run checks on N
    if asserts:
        assert N.min() >= 0
        assert True not in np.isnan(N)
        assert True not in np.isinf(N)
    #add debug plot
    if debug:
        dp._next_plot('N', N)

    #get K raster
    K = pop.K
    #run checks on K
    if asserts:
        assert K.min() >= 0
        assert True not in np.isnan(K)
        assert True not in np.isinf(K)
    #add debug plot
    if debug:
        dp._next_plot('K', K)
    #NOTE: SMALL VALUES ON THE K RASTER TOTALLY SCREW UP DENSITY-DEPENDENCE
    #BECAUSE THEY GENERATE HUGE NEGATIVE dN/dt VALUES BY DIVISION; BRAINSTORM
    #SOME MORE SOPHISTICATED, REALISTIC WAY TO HANDLE!
    #Perhaps there can be a separate step after calculating dNdt where any
    #0-cells are converted from the extremely negative numbers calcuated
    #there to simply the negative of the number of individuals currently found
    #there (n) minus n*cell value, so that the probability of individuals
    #there dying winds up approaching 1?... something like that

    #calc dNdt
    dNdt = _calc_dNdt(land = land, R = pop.R, N = N, K = K,
                                    pop_growth_eq = 'logistic')
    #run checks on dNdt
    if asserts:
        assert True not in np.isnan(dNdt)
        assert True not in np.isinf(dNdt), ('The following cells are '
            'infinite: \n\t%s') % str([i for i, n in enumerate(
                                        dNdt.ravel()) if np.isinf(n)])
    #add debug plot
    if debug:
        dp._next_plot('dNdt', dNdt)

    #calculate N_b (raster of estimated births)
    N_b = _calc_N_b(b = pop.b,
        n_births_distr_lambda = pop.n_births_distr_lambda, n_pairs = n_pairs)
    #run checks on N_b
    if asserts:
        assert N_b.min() >= 0
        assert True not in np.isnan(N_b)
        assert True not in np.isinf(N_b)
    #add debug plot
    if debug:
        dp._next_plot('N_b', N_b)

    #calc N_d (raster of deaths)
    N_d = _calc_Nd(N_b = N_b, dNdt = dNdt)
    #run checks on N_d
    if asserts:
        assert True not in np.isnan(N_d)
        assert True not in np.isinf(N_d)
    #add debug plot
    if debug:
        dp._next_plot('N_d', N_d)

    #calc d (raster of probabilities of density-dependent death)
    d = _calc_d(N_d = N_d, N = N, d_min = pop.d_min, d_max = pop.d_max)
    #run checks on d
    if asserts:
        assert d.min() >= 0, 'd.min() is %0.2f, at %s' % (d.min(),
                                                        str(d.argmin()))
        assert d.max() <= 1, 'd.max() is %0.2f' % d.max()
        assert True not in np.isnan(d)
        assert True not in np.isinf(d)
    #add debug plot
    if debug:
        dp._next_plot('d', d)

    #If births should happen after (and thus not be included in the 
    #calculation of) deaths, then instead of having mated and dispersed
    #babies up above, do it now (i.e. now that the d raster has been 
    #calculated)
    #(see following 2 NOTE)
    #NOTE: 04/28/17: LEAVING THE FOLLOWING NOTE FOR NOW, BUT I ADDED
    #THIS EARLIER, BEFORE TODAY ADDING THE births_before_deaths ARGUMENT 
    #THAT DEFAULTS TO False. The problem I overlooked in the steps that I
    #laid out in the following NOTE is that I failed to realize that the
    #probability of a baby's death
    #would be based on the density of parents before it was born, such
    #that deaths of babies would actually be artificially low, and I believe
    #this is what was leading to my constant, inevitable exponential population
    #dynamics
    #NOTE: The way I'm looking at it now, calculating likely births from number
    #of existing pairs, then calculating death raster, then having births take
    #place, then having per-individual death probabilies determined and having
    #deaths take place (such that newborns could also immediately die after
    #birth, based on the death raster calculated before they were born, and
    #even based on individs who weren't their parents, assuming they dispersed
    #across cell lines), which is what I'm doing currently, is the only
    #obvious/reasonable way to do this. Of course there will always be an
    #inherent rub between reality and discrete-time simulations because of
    #the continuous/discrete time conflict. But in this case, is there any
    #strong reason to CHANGE THIS? Or are there any foreseeable and undesirable
    #results/effects?... NEED TO PUT MORE THOUGHT INTO THIS LATER.
    if not births_before_deaths:
        #Feed the land and mating pairs to the mating functions, to produce
        #and disperse zygotes
        pop._do_mating(pairs, burn)
    #Get death probabilities
    death_probs = d[pop.cells[:,1], pop.cells[:,0]]
    #If with_selection (i.e. if death probs should account for fitness),
    #then use the d raster and individuals' fitnesses to calculate
    #per-individual probabilities of death
    if with_selection:
        death_probs = selection._calc_prob_death(pop, death_probs)
    #run checks on death_probs
    if asserts:
        assert np.alltrue(death_probs >= 0)
        assert np.alltrue(death_probs <= 1)

    #kill (and track kills) based on max_age
    num_killed_age = 0
    if pop.max_age is not None:
        death_probs[pop._get_age() > pop.max_age] = 1
        num_killed_age = np.sum(death_probs == 1)

    #Use the per-individual death probabilities to carry out mortality 
    num_deaths = _do_mortality(land, pop, death_probs)
    pop._set_coords_and_cells()
    pop.n_deaths.append(num_deaths)

    #Check if extinct, and return result
    extinct = pop._check_extinct()
    return(extinct)


