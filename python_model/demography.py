#!/usr/bin/python
#demography.py

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy import random as r
from scipy import interpolate
from scipy.spatial import cKDTree
from sklearn.preprocessing import normalize
from collections import OrderedDict as OD

import landscape
import mating
import selection


'''Functions to control demography/population dynamics.'''


#TODO:
    # - probably soon get rid of all the debug and print statements in the pop_dynamics function soon!
    # - also probably get rid of any optional arguments fed into if statements that I wound up not using?


def logistic_eqxn(R, N, K):
    dNdt = R*(1-(N/K))*N 
    return(dNdt)


def logistic_soln(x0, R, t):
        return(1/(1+((1/x0)-1)*np.e**(-1*R*t)))



def calc_dNdt(land, N, K, params, pop_growth_eq = 'logistic'):

    dims = land.dims
    
    #if using the logistic equation for pop growth
    #NOTE: For now this is the only option and the default, but could easily offer other equations later if desired
    if pop_growth_eq == 'logistic':
        #use logistic eqxn, with pop intrinsic growth rate from params['r'], to generate current growth-rate raster
        dNdt = logistic_eqxn(params['R'], N, K)
    
    return(landscape.Landscape(dims, dNdt))




def kill(land, pop, params, death_probs):
    deaths = np.array(list(death_probs.keys()))[np.bool8(r.binomial(n = 1, p = list(death_probs.values())))]
    [land.mating_grid.remove(pop.individs[ind]) for ind in deaths]
    [pop.individs.pop(ind) for ind in deaths]

    return(len(deaths))




def pop_dynamics(land, pop, params, with_selection = True, burn = False, age_stage_d = None, births_before_deaths = False, debug = None):
    '''Generalized function for implementation population dynamics. Will carry out one round of mating and
    death, according to parameterization laid out in params dict.


       If with_selection == False, only density-dependent death will occur.

       If age_stage_d != None, an nx2 Numpy array should be provided, with column 1 listing the ages or
       stages and column 2 stipulating the relative probabilites of death for each age or stage (i.e. must be
       a discrete probability distribution, i.e. sum to 1). NOTE: THIS ISN'T IMPLEMENTED YET!


       If burn_in == True, selection will be coerced to False, the Population.K attribute will be updated,
       stationarity metrics will be assessed, and the function will return the decision it has reached
       regarding stationarity. NOTE: THIS ISN'T IMPLEMENTED YET!

    '''

    ##########################################
    ####DEBUG OPTION, WHICH WILL PLOT THE RASTERS FOR EACH STEP OF THE PROCESS:
    #debug = False
    #debug = True
    ##########################################



    #NOTE: Would be good to write in option to ONLY mate or ONLY die when the function is run (so that a pop could be run with mating every x timesteps but mortality each timestep, for example)



    ######calc num_pairs raster (use the calc_pop_density function on the centroids of the mating pairs)
    pairs = pop.find_mating_pairs(land, params)


    p_x = [float(np.mean((pop.individs[pairs[i,0]].x, pop.individs[pairs[i,1]].x))) for i in range(pairs.shape[0])]
    p_y = [float(np.mean((pop.individs[pairs[i,0]].y, pop.individs[pairs[i,1]].y))) for i in range(pairs.shape[0])]

    #use the Landscape_Stack.density_grid_stack to calculate a raster of the expected number of pairs
    #NOTE: because the window_width of that object is much wider than mu_dispersal in the parameterizations
    #I've typically been running, this produces much wider areas of high values in the 
    #N_b, N_d, and hence d rasters than the old method did (which had its window_width set to mu_dispersal,
    #leading of course to huge lags in the old calc_density algorithm)
    #need to think about whether I need to attend to this, and if so, how
    n_pairs = np.clip(land.density_grid_stack.calc_density(p_x, p_y), a_min = 0, a_max = None)

    n_pairs[np.isnan(n_pairs)] = 0
    assert n_pairs.min() >= 0, 'n_pairs.min() == %0.2f' %(n_pairs.min())  

    #NOTE: I anticipate that using mu_dispersal as the density-calculating window should produce a slightly more realistic expectation
    #of the number of births per cell in a few steps; using max(1, mu_dispersal) to avoid
    #window_lengths <1 (needlessly computationally expensive)

    assert True not in np.isnan(n_pairs)
    assert True not in np.isinf(n_pairs)

    if debug != None:
        fig = plt.figure()
        plt.suptitle('TIMESTEP: %i' % debug)
        var = 'n_pairs'
        plt.subplot(241)
        plt.title(var)
        pop.show(land,1, markersize = 8)
        plt.imshow(eval(var), interpolation = 'nearest', cmap = 'terrain')
        plt.colorbar()
        print(var)
        print(eval(var))
        #raw_input()




    ######if births should happen before (and thus be included in the calculation of) deaths, then mate and disperse babies now
    if births_before_deaths == True:
        #Feed the mating pairs and params['b'] to the mating functions, to produce and disperse zygotes
        pop.mate(land, params, pairs, burn)




    ######calc N raster, set it as pop.N, then get it  
    #NOTE: 04/28/17 this way the N used to calculate density-dependent deaths this turn incorporates the babies born this turn, if births_before_deaths == True
    pop.calc_density(land, set_N = True)
    N = pop.N.raster
    assert N.min() >= 0
    assert True not in np.isnan(N)
    assert True not in np.isinf(N)


    if debug != None:
        var = 'N'
        plt.subplot(242)
        plt.title(var)
        pop.show(land,1, markersize = 8)
        plt.imshow(eval(var), interpolation = 'nearest', cmap = 'terrain')
        plt.colorbar()
        print(var)
        print(eval(var))
        #raw_input()



    ######get K raster
    K = pop.K.raster
    #K[K==0]=1  #values of 0 produce np.inf by division, so must avoid this
    #NOTE: Plotted dNdt against N for different values of K, to see how strong density-dependent
    #selection would be for individs landing in lowest-K cells. 1 seemed still very strong, but a bit more
    #permissive than the absurdly negative values that were produced by adding, say, 10e-1 or less; for now
    #leaving 1 as an arbitrary lower bound on K

    #if params dict has a K_cap parameter, use this to set the max cell value for the K raster
    #if 'K_cap' in params.keys():
        #K[K>float(params['K_cap'])] = float(params['K_cap'])
        #print(K.max())
    assert K.min() >= 0
    assert True not in np.isnan(K)
    assert True not in np.isinf(K)


    if debug != None:
        var = 'K'
        plt.subplot(243)
        plt.title(var)
        pop.show(land,1, markersize = 8)
        plt.imshow(eval(var), interpolation = 'nearest', cmap = 'terrain')
        plt.colorbar()
        print(var)
        print(eval(var))
        #raw_input()



    #NOTE: SMALL VALUES ON THE K RASTER TOTALLY SCREW UP DENSITY-DEPENDENCE BECAUSE 
    #THEY GENERATE HUGE NEGATIVE dN/dt VALUES BY DIVISION; BRAINSTORM SOME MORE SOPHISTICATED, REALISTIC WAY TO HANDLE!
    #Perhaps there can be a separate step after calculating dNdt where any  0-cells
    #are convereted from the extremely negative numbers calcuated there to simply the
    #negative of the number of individuals currently found there (n) minus n*cell
    #value, so that the probability of individuals there dying winds up approaching
    #1?... something like that



    ######use N, K, and params['r'] to calc dN/dt raster, according to the chosen pop growth eqxn
        #NOTE: Currently, only option and default is basic logistic eqxn: dN/dt = r*(1-N/K)*N
    with np.errstate(divide='ignore', invalid='ignore'):
        dNdt = calc_dNdt(land, N, K, params, pop_growth_eq = 'logistic').raster

    #NOTE: The next line used to replace each cell in dNdt where K<1 with that cell's -N value. But it was
    #really slow to run, and at those cells where K ~= 0 it really won't make a difference how negative the
    #dNdt values there are, so the following line just makes this much simpler.
    dNdt = np.clip(dNdt, a_min = -1*N.max(), a_max = None)

    assert True not in np.isnan(dNdt)
    assert True not in np.isinf(dNdt), 'The following cells are infinite: \n\t%s' % str([i for i, n in enumerate(dNdt.ravel()) if np.isinf(n)])

    if debug != None:
        var = 'dNdt'
        plt.subplot(244)
        plt.title(var)
        pop.show(land,1, markersize = 8)
        plt.imshow(eval(var), interpolation = 'nearest', cmap = 'terrain')
        plt.colorbar()
        print(var)
        print(eval(var))
        print(eval(var).min())
        print(eval(var).argmin())
        #raw_input()



    ######use n_pairs, params['b'] (i.e. probability that a pair gives birth), and params['lambda_offspring'] (i.e.
    #expected number of births per pair) to calculate the N_b raster (expected number of births in each cell)
    #according to N_b = b*E[N_b/pair]*n_pairs where E[N_b/pair] = lambda (b/c births Poission distributed)
        #NOTE: Would be good to in the future to use neg_binom dist to allow births to be modeled as Poisson
        #(using correct parameterization of neg_binom to approximate that) or overdispersed (to be parmeterized with field data)

        #NOTE: Does it matter that I'm not attempting to account for the fact that pairs will disperse out of
        #the cell? Almost certainly a larger effect with increasing dispersal-kernel distance. But nonetheless,
        #would it even make sense (and even BE TRACTABLE) to calculate it any other way? My intuition is NO, but 
        #worth considering further later on...

        #NOTE: Currently, birth rate is intrinsic to the species and cannot be made a function of density or
        #spatial fitness. Might be worth adding this at some point in the future though.
    N_b = params['b'] * params['lambda_offspring'] * n_pairs
    assert N_b.min() >= 0
    assert True not in np.isnan(N_b)
    assert True not in np.isinf(N_b)
 
    if debug != None:
        var = 'N_b'
        plt.subplot(245)
        plt.title(var)
        pop.show(land,1, markersize = 8)
        plt.imshow(eval(var), interpolation = 'nearest', cmap = 'terrain')
        plt.colorbar()
        print(var)
        print(eval(var))
        #raw_input()

  

    ######Use N_b - dN/dt to calculate N_d raster (expected number of deaths in each cell)
    N_d = N_b - dNdt
    assert True not in np.isnan(N_d)
    assert True not in np.isinf(N_d)

    if debug != None:
        var = 'N_d'
        plt.subplot(246)
        plt.title(var)
        pop.show(land,1, markersize = 8)
        plt.imshow(eval(var), interpolation = 'nearest', cmap = 'terrain')
        plt.colorbar()
        print(var)
        print(eval(var))
        #raw_input()

   

    ######Use N_d/N to calculate d raster (probability of death of individuals in each cell)
    #d = N_d/N
    #d[d<0] = 0

    #Instead of the above, just normalize N_d to use it as the raster of the probability of death of individuals in each cell
    #NOTE: THIS ACTUALLY SEEMS LIKE A BAD APPROACH, BECAUSE IT FORCES THE SPATIAL DISTRIBUTION OF d TO MIRROR
    #THE SPATIAL DISTRIBUTION OF RELATIVE dNdt AND N_b VALUES, but in reality d should be similarly high in low-K
    #cells with low N and hi-K cells with hi N
    #d = (N_d - N_d.min())/(N_d.max() - N_d.min())
    #NOTE: Instead, calculate d as simply N_d/N, then control for some numerical artefacts that an inappropriate

    #floor N_d to 0, since negative N_d values are not going to generate new individuals at this point, but will create wonky values in d
    #N_d[N_d <= 0] = 0.0000001  
    #divide N_d/N
    d = N_d/N
    #fix infinties if they arise (where N == 0)
    d[np.isinf(d)] = 0
    d[np.isnan(d)] = 0
    d[d<0] = 0
    #constrain to the min and max d values
    d_min = params['d_min']
    d_max = params['d_max']
    d[d<d_min] = d_min
    d[d>d_max] = d_max


    #d = normalize(N_d, norm = 'l2') 
    #NOTE: the sklearn normalize function seems to tweak the distribution of relative
    #values a bit (but JUST a bit), so I should read about the normalization methods in detail later, but for
    #now using this instead of the above approach because the above approach is guaranteed to produce both 0
    #and 1 values, and particularly the 1s are problematic
    #NOTE: For some reason I was still winding up with slightly negative values after running normalize...? So
    #just coerce them to 0
    #d[d<0] = 0
    
    #Actually, in fact just coercing everything below a min to that value, and everything above a max to
    #that value, just to impose an arbitrary min and max mortality probability. (Feels wrong to have
    #certainty of death or survival in some locations.) So in that case, stick with my max-diff normalization,
    #then coercing the ends to the min and max probs
    #d[d<d_min] = d_min
    #d[d>d_max] = d_max


    assert d.min() >= 0, 'd.min() is %0.2f, at %s' % (d.min(), str(d.argmin()))
    assert d.max() <= 1, 'd.max() is %0.2f' % d.max()
    assert True not in np.isnan(d)
    assert True not in np.isinf(d)

    if debug != None:
        var = 'd'
        plt.subplot(247)
        plt.title(var)
        pop.show(land,1, markersize = 8)
        plt.imshow(eval(var), interpolation = 'nearest', cmap = 'terrain')
        plt.colorbar()
        print(var)
        print(eval(var))
        raw_input('\n\n\nPress <Enter> for next timestep...\n\n')




    ######If births should happen after (and thus not be included in the calculation of) deaths, then instead
         #of having mated and dispersed babies up above, do it now (i.e. now that the d raster has been calculated)
         #(see following 2 NOTE)
        
        #NOTE: 04/28/17: LEAVING THE FOLLOWING NOTE FOR NOW, BUT I ADDED THIS EARLIER, BEFORE TODAY ADDING THE
        #births_before_deaths ARGUMENT THAT DEFAULTS TO False. The problem I overlooked in the setps that I
        #laid out in the following NOTE is that I failed to realize that the probability of a baby's death
        #would be based on the density of parents before it was born, such that deaths of babies would actually
        #be artificially low, and I believe this is what was leading to my constant, inevitable exponential
        #population dynamics
            #NOTE: The way I'm looking at it now, calculating likely births from number of existing pairs, then
            #calculating death raster, then having births take place, then having per-individual death probabilies
            #determined and having deaths take place (such that newborns could also 
            #immediately die after birth, based on the death raster calculated before they were born, and even
            #based on individs who weren't their parents, assuming they dispersed across cell lines), which is what
            #I'm doing currently, is the only obvious/reasonable way to do this. Of course there will always be an
            #inherent rub between reality and discrete-time simulations because of the continuous/discrete time
            #conflict. But in this case, is there any strong reason to CHANGE THIS? Or are there any foreseeable
            #and undesirable results/effects?... NEED TO PUT MORE THOUGHT INTO THIS LATER.

    if births_before_deaths == False:
        #Feed the mating pairs and params['b'] to the mating functions, to produce and disperse zygotes
        pop.mate(land, params, pairs, burn)



    ######Now implement deaths

    #If with_selection = True, then use the d raster and individuals' relative fitnesses to calculate
    #per-individual probabilities of death
    #NOTE: FOR NOW JUST USING A SINGLE LOCUS, BUT LATER NEED TO IMPLEMENT FOR AN ARBITRARY NUMBER OF LOCI, AND A MAP
    #OF THOSE LOCI ONTO TRAITS
    if with_selection == True:
    
        death_probs = selection.get_prob_death(pop, {i:d[int(ind.y), int(ind.x)] for i, ind in pop.individs.items()})

        #print(death_probs.items()[0:20])
        #print('VS')
        #print(dict([(i, d[int(ind.y), int(ind.x)]) for i,ind in pop.individs.items()]).items()[0:20])

    elif with_selection == False:
        death_probs = OD({i:d[int(ind.y), int(ind.x)] for i, ind in pop.individs.items()})
        assert np.alltrue(np.array(list(death_probs.values())) >= 0)
        assert np.alltrue(np.array(list(death_probs.values())) <= 1)

    
    if params['island_val'] > 0:
        death_probs.update({i:1 for i,v in pop.get_habitat_by_land_ind(scape_num = land.n_island_mask_scape).items() if v})
        num_killed_isle = len({i:1 for i,v in pop.get_habitat_by_land_ind(scape_num = land.n_island_mask_scape).items() if v})
        print('\n\tINDIVIDS KILLED OUTSIDE ISLANDS: %i  (%0.3f%% of pop)\n' % (num_killed_isle, num_killed_isle/pop.Nt[::-1][0]))
        
    

    #If age_stage_d != None then use the age_stage_d array, and pop.get_age(), and the d raster 
    #to calculate per-individual death probabilities
    if age_stage_d != None:
        pass
        

    #Feed the per-individual death probabilities into the kill function, which will probabilistically generate deaths
    #and cull those individuals, and will return the number of deaths
    num_deaths = kill(land, pop, params, death_probs)
    pop.n_deaths.append(num_deaths)

    print('\n\t%i individuals dead' % num_deaths)



    #Check if extinct
    extinct = pop.check_extinct()
        


    return(extinct)
    #If burn_in == True:
    #if burn_in == True:
    #    pass
        #Update pop.K
        #Run the burn_in functions to update the running metrics, and assess stationarity, and return decision
        #of whether or not stationarity has been achieved

        #return decision
        
    #Else exit
    #else:
        #pass




