#!/usr/bin/python
#demography.py

'''
Functions to implement demographic operations (birth and death).
'''

#geonomics imports
from geonomics.ops.selection import _calc_prob_death
from geonomics.utils.viz import _check_display

#other imports
import numpy as np
import matplotlib as mpl
_check_display()
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
    def __init__(self, land, spp, timestep):
        self.subplot = 241
        self.timestep = timestep
        self.spp = spp
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
            self.spp._plot(land = self.land, lyr_num = 0, size = 8, hide_land = True)
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

def _calc_n_pairs(pairs, spp):
    #if there are pairs
    if pairs.size > 0:
        #get their coordinates
        pairs_coords = spp._get_coords(individs = pairs.flatten())
        #take just the x and y coordinates from pairs_coords, reshape them to
        #match pairs.shape (i.e.nx2, where n = len(spp)), and then mean across
        #axis 1 (to get the pair's midpoint
        #x and y coordinates)
        p_x = pairs_coords[:,0].reshape(pairs.shape).mean(axis = 1)
        p_y = pairs_coords[:,1].reshape(pairs.shape).mean(axis = 1)

        #use the spp._dens_grids to calculate a raster of the expected number
        #of pairs
        #NOTE: because the window_width of that object is much wider than 
        #dispersal_mu in the parameterizations I've typically been running, 
        #this produces much wider areas of high values in the N_b, N_d, and
        #hence d rasters than the old method did (which had its window_width
        #set to dispersal_mu, leading of course to huge lags in the old 
        #calc_density algorithm) need to think about whether I need to attend
        #to this, and if so, how
        n_pairs = np.clip(spp._dens_grids._calc_density(p_x, p_y), a_min = 0,
                                                                a_max = None)

        #remove NaNs
        n_pairs[np.isnan(n_pairs)] = 0

    #otherwise just return a 0 raster
    else:
        n_pairs = np.zeros(spp._land_dim)

    return n_pairs


#the logistic eqxn
def _calc_logistic_growth(R, N, K):
    dNdt = R*(1-(N/K))*N
    return dNdt


def _calc_logistic_soln(x0, R, t):
        return 1/(1+((1/x0)-1)*np.e**(-1*R*t))


def _calc_dNdt(R, N, K, pop_growth_eq = 'logistic'):
    with np.errstate(divide='ignore', invalid='ignore'):
        #if using the logistic equation for pop growth
        #NOTE: For now this is the only option and the default, but could 
        #easily offer other equations later if desired
        if pop_growth_eq == 'logistic':
            #use logistic eqxn, with spp intrinsic growth rate (spp.R)
            #to generate current growth-rate raster
            dNdt = _calc_logistic_growth(R, N, K)
    #coerce extreme negative values, NaNs, and infs to a reasonable value
    #(extreme negative values are produced in areas where the K raster
    #has extremely small values)
    dNdt = np.clip(dNdt, a_min = -1*N.max(), a_max = None)
    dNdt[np.isnan(dNdt)] = -1*N.max()
    dNdt[np.isinf(dNdt)] = -1*N.max()
    return dNdt


def _calc_N_b(b, n_births_distr_lambda, n_pairs):
    #Use n_pairs, spp.b (i.e. probability that a pair gives birth),
    #and spp.n_births_distr_lambda (i.e. expected number of births
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
    # NOTE: ignore x/0. and 0./0. warnings
    with np.errstate(divide='ignore', invalid='ignore'):
        d = N_d/N
    #fix infinties and NaNs and negatives if they arise
    #(they occur where N ==0 or where N_d == 0)
    d[np.isnan(d)] = 0 
    #constrain to the min and max d values
    d = np.clip(a = d, a_min = d_min, a_max = d_max)

    #NOTE: DEH: 12-08-18: Got rid of the infinity check, because it was
    #coercing all infinities to 0, but in many cases infinities should actually
    #have been 1s, so it created numerous problems that just clipping to d_min
    #and d_max takes care of anyhow; commented out the NaN check below, but not
    #entirely clear that I don't need it...
    #d[np.isnan(d)] = 0
    return d


def _do_mortality(spp, death_probs):
    deaths = np.array([*spp])[np.bool8(r.binomial(n = 1, p = death_probs))]
    if len(deaths) > 0:
        ig = itemgetter(*deaths)
        [spp.pop(ind) for ind in deaths];
    return len(deaths)


def _do_pop_dynamics(spp, land, with_selection = True, burn = False,
    births_before_deaths = False, asserts = True, debug = False):
    '''Generalized function for implementation population dynamics.
    Will carry out one round of mating and death, according to parameterization
    laid out in params dict (which were grabbed as Species attributes).

       If with_selection == False, only density-dependent death will occur.

       If burn == True, selection will be coerced to False, the Species.K
       attribute will be updated, stationarity metrics will be assessed, and
       the function will return the decision it has reached regarding
       stationarity.

    '''

    #NOTE: Would be good to write in option to ONLY mate or ONLY die when the
    #function is run (so that a spp could be run with mating every x timesteps
    #but mortality each timestep, for example)

    #create the debug-plotter, id debug == True
    if debug:
        timestep = -99999
        dp = _DebugPlotter(land, spp, timestep)

    #find mating pairs
    pairs = spp._find_mating_pairs()

    #calc num_pairs raster (use the calc_pop_density function on the
    #midpoints of the mating pairs)
    n_pairs = _calc_n_pairs(pairs = pairs, spp = spp)
    #run checks on n_pairs
    if asserts:
        assert n_pairs.min() >= 0, 'n_pairs.min() == %0.2f' %(n_pairs.min())
        assert not np.any(np.isnan(n_pairs))
        assert not np.any(np.isinf(n_pairs))
    #add debug plot
    if debug:
        dp._next_plot('n_pairs', n_pairs)

    #if births should happen before (and thus be included in the calculation of)
    #deaths, then mate and disperse babies now
    if births_before_deaths:
        #Feed the land and mating pairs to spp.do_mating, to produce and
        #disperse zygotes
        spp._do_mating(land, pairs, burn)

    #calc N raster, set it as spp.N, then get it  
    spp._calc_density(set_N = True)
    N = spp.N
    #run checks on N
    if asserts:
        assert N.min() >= 0
        assert not np.any(np.isnan(N))
        assert not np.any(np.isinf(N))
    #add debug plot
    if debug:
        dp._next_plot('N', N)

    #get K raster
    K = spp.K
    #run checks on K
    if asserts:
        assert K.min() >= 0
        assert not np.any(np.isnan(K))
        assert not np.any(np.isinf(K))
    #add debug plot
    if debug:
        dp._next_plot('K', K)

    #calc dNdt
    dNdt = _calc_dNdt(R = spp.R, N = N, K = K,
                                    pop_growth_eq = 'logistic')
    #run checks on dNdt
    if asserts:
        assert not np.any(np.isnan(dNdt))
        assert not np.any(np.isinf(dNdt)), ('The following cells are '
            'infinite: \n\t%s') % str([i for i, n in enumerate(
                                        dNdt.ravel()) if np.isinf(n)])
    #add debug plot
    if debug:
        dp._next_plot('dNdt', dNdt)

    #calculate N_b (raster of estimated births)
    N_b = _calc_N_b(b = spp.b,
        n_births_distr_lambda = spp.n_births_distr_lambda, n_pairs = n_pairs)
    #run checks on N_b
    if asserts:
        assert N_b.min() >= 0
        assert not np.any(np.isnan(N_b))
        assert not np.any(np.isinf(N_b))
    #add debug plot
    if debug:
        dp._next_plot('N_b', N_b)

    #calc N_d (raster of deaths)
    N_d = _calc_Nd(N_b = N_b, dNdt = dNdt)
    #run checks on N_d
    if asserts:
        assert not np.any(np.isnan(N_d))
        assert not np.any(np.isinf(N_d))
    #add debug plot
    if debug:
        dp._next_plot('N_d', N_d)

    #calc d (raster of probabilities of density-dependent death)
    d = _calc_d(N_d = N_d, N = N, d_min = spp.d_min, d_max = spp.d_max)

    #run checks on d
    if asserts:
        assert d.min() >= 0, 'd.min() is %0.2f, at %s' % (d.min(),
                                                        str(d.argmin()))
        assert d.max() <= 1, 'd.max() is %0.2f' % d.max()
        assert not np.any(np.isnan(d))
        assert not np.any(np.isinf(d))
    #add debug plot
    if debug:
        dp._next_plot('d', d)

    #If births should happen after (and thus not be included in the 
    #calculation of) deaths, then instead of having mated and dispersed
    #babies up above, do it now (i.e. now that the d raster has been 
    #calculated)
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
        spp._do_mating(land, pairs, burn)
    #Get death probabilities
    death_probs = d[spp.cells[:,1], spp.cells[:,0]]
    #If with_selection (i.e. if death probs should account for fitness),
    #then use the d raster and individuals' fitnesses to calculate
    #per-individual probabilities of death
    if with_selection:
        death_probs = _calc_prob_death(spp, death_probs)
    #run checks on death_probs
    if asserts:
        assert np.alltrue(death_probs >= 0)
        assert np.alltrue(death_probs <= 1)

    #kill (and track kills) based on max_age
    num_killed_age = 0
    if spp.max_age is not None:
        death_probs[spp._get_age() > spp.max_age] = 1
        num_killed_age = np.sum(death_probs == 1)

    #Use the per-individual death probabilities to carry out mortality 
    num_deaths = _do_mortality(spp, death_probs)
    spp._set_coords_and_cells()
    spp.n_deaths.append(num_deaths)

    #Check if extinct, and return result
    extinct = spp._check_extinct()
    return(extinct)


