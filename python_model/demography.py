#!/usr/bin/python
#demography.py

import numpy as np
import matplotlib.pyplot as plt
from numpy import random as r
from scipy import interpolate
from scipy.spatial import cKDTree
from sklearn.preprocessing import normalize

import landscape
import mating


'''Functions to control demography/population dynamics.'''


def calc_pop_density(land, coords, window_width = None, normalize_by = 'none', min_0 = True, max_1 = False, max_val = None):

    '''
    Calculate an interpolated raster of local population density, using a window size of window_width.
    Valid values for normalize_by currently include 'census' and 'none'. If normalize_by = 'census', max_1 =
    True will cause the output density raster to vary between 0 and 1, rather than between 0 and the current
    max normalized density value. Window width will default to 1/10 of the larger raster dimension.
    '''

    #window width default to 1/10 the maximum landscape dimension
    if window_width == None:
        window_width = max(land.dims)*0.1

    #shorthand
    dims = land.dims

    #get a list of pop's coord-tuples
    c = coords

    #make window_width a float, to avoid Py2 integer-division issues
    window_width = float(window_width)
    
    #create meshgrid using window_width/2 as step size
    grid_j, grid_i = np.mgrid[0:dims[0]:complex("%ij" % (dims[0]/(window_width/2))), 0:dims[1]:complex("%ij" % (dims[1]/(window_width/2)))]

    #flatten the arrays, so that I can run over them in a single for loop
    gj = grid_j.ravel()
    gi = grid_i.ravel()

    #make lists of tuples, of same length as gj, containing the window ll and ur coords
    window_ll = [(max(gj[n]-window_width/2, 0), max(gi[n]-window_width/2, 0)) for n in range(len(gj))]   #constrain min window vals to 0
    window_ur = [(min(gj[n]+window_width/2, land.dims[0]), min(gi[n]+window_width/2, land.dims[1])) for n in range(len(gj))] #constrain max window vals to each respective land dimension
    assert len(window_ll) == len(gj)
    assert len(window_ur) == len(gj)

    #make a list of the sizes of each window
    window_size = [(window_ur[n][0] - window_ll[n][0]) * (window_ur[n][1] - window_ll[n][1]) for n in range(len(gj))]#calculate size of this window (not always the same because of edge effects
    assert len(window_size) == len(gj)
    
    #make a list of the counts of all individs within each window
    window_ct = [len([ind for ind in range(len(c)) if (c[ind][0]>window_ll[n][0] and c[ind][0]<=window_ur[n][0]) and (c[ind][1]>window_ll[n][1] and c[ind][1]<=window_ur[n][1])]) for n in range(len(gj))] 
    assert len(window_ct) == len(gj)

    #divide window counts by window sizes
    window_dens = [window_ct[n]/window_size[n] for n in range(len(window_ct))] #divide by window size
    assert len(window_dens) == len(gj)

    #if normalize_by == census, then divide each density by total pop census size
    if normalize_by == 'census':
        N = pop.census()
        window_dens = [dens/N for dens in window_dens]
    elif normalize_by == 'none':
        pass

    else:  #POTENTIALLY ADD OTHER OPTIONS HERE, i.e. to normalize by starting size?
        pass 

    #interpolate resulting density vals to a grid equal in size to the landscape
    new_gj, new_gi = np.mgrid[0:dims[0]:complex("%ij" % (dims[0])), 0:dims[1]:complex("%ij" % (dims[1]))]
    dens = interpolate.griddata(zip(list(gj), list(gi)), window_dens, (new_gj, new_gi), method = 'cubic')



    if normalize_by <> 'none':

        #if max_1 == True, set max_val to dens.max(), such that the density raster output will be normalized to
        #its own max, and thus vary between 0 and 1; else set to 1, and the output raster will vary between 0 and the current max value
        if max_1 == True:
            max_val = dens.max()
        elif max_1 == False:
            max_val = 1

        #Use max_val to normalize the density raster to either 0 to its current max val or
        #0 to 1, to make sure the interpolation didn't generate any values slightly outside this range
        norm_factor = max_val - dens.min()
        dens = (dens - dens.min())/norm_factor

    else:
        pass

    if min_0 == True:
        dens[dens<0] = 0

    if max_val <> None:
        dens[dens>max_val] = max_val
    
    return(landscape.Landscape(dims, dens))





def logistic_eqxn(r, N, K):
    dNdt = r*(1-(N/K))*N
    return(dNdt)



def calc_dNdt(land, N, K, params, pop_growth_eq = 'logistic'):

    dims = land.dims
    
    #if using the logistic equation for pop growth
    #NOTE: For now this is the only option and the default, but could easily offer other equations later if desired
    if pop_growth_eq == 'logistic':
        #use logistic eqxn, with pop intrinsic growth rate from params['r'], to generate current growth-rate raster
        dNdt = logistic_eqxn(params['r'], N, K)
    
    return(landscape.Landscape(dims, dNdt))




def mort(land, pop, params, death_probs):
    deaths = [i for i, p in death_probs.items() if bool(r.binomial(1, p, 1)) == True]

    [pop.individs.pop(ind) for ind in deaths]

    return(len(deaths))




def pop_dynamics(land, pop, params, selection = True, burn_in = False, age_stage_d = None, d_min = 0.01, d_max = 0.99):
    '''Generalized function for implementation population dynamics. Will carry out one round of mating and
    deaths, according to parameterization laid out in params dict.


       If selection == False, only density-dependent death will occur.

       If age_stage_d <> None, an nx2 Numpy array should be provided, with column 1 listing the ages or
       stages and column 2 stipulating the relative probabilites of death for each age or stage (i.e. must be
       a discrete probability distribution, i.e. sum to 1). NOTE: THIS ISN'T IMPLEMENTED YET!


       If burn_in == True, selection will be coerced to False, the Population.K attribute will be updated,
       stationarity metrics will be assessed, and the function will return the decision it has reached
       regarding stationarity. NOTE: THIS ISN'T IMPLEMENTED YET!

    '''

    #NOTE: Would be good to write in option to ONLY mate or ONLY die when the function is run (so that a pop could be run with mating every x timesteps but mortality each timestep, for example)



    #calc num_pairs raster (use the calc_pop_density function on the centroids of the mating pairs)
    pairs = pop.find_mating_pairs(pop, params)
    print(pairs)
    p_x = [int(np.mean((pop.individs[pairs[i,0]].x, pop.individs[pairs[i,1]].x))) for i in range(pairs.shape[0])]
    p_y = [int(np.mean((pop.individs[pairs[i,0]].y, pop.individs[pairs[i,1]].y))) for i in range(pairs.shape[0])]
    n_pairs = calc_pop_density(land, zip(p_x, p_y), max(1, params['mu_dispersal']), min_0 = True).raster
    n_pairs[np.isnan(n_pairs)] = 0
    plt.imshow(n_pairs)
    plt.colorbar()
    assert n_pairs.min() >= 0, 'n_pairs.min() == %0.2f' %(n_pairs.min())
    assert True not in np.isnan(n_pairs)
    assert True not in np.isinf(n_pairs)


    #NOTE: I anticipate that using mu_dispersal as the density-calculating window should produce a slightly more realistic expectation
    #of the expected number of births per cell in a few steps; using max(1, mu_dispersal) because small
    #window_lengths are needlessly computationally expensive

    #get N raster
    N = pop.N.raster
    assert N.min() >= 0
    assert True not in np.isnan(N)
    assert True not in np.isinf(N)


    #print('Plotting N')
    #plt_N = landscape.Landscape(land.dims, N)
    #pop.show(plt_N)
    #raw_input()


    #get K raster
    K = pop.K.raster
    K[K==0]=1  #values of 0 produce np.inf by division, so must avoid this
    #NOTE: Plotted dNdt against N for different values of K, to see how strong density-dependent
    #selection would be for individs landing in lowest-K cells. 1 seemed still very strong, but a bit more
    #permissive than the absurdly negative values that were produced by adding, say, 10e-1 or less; for now
    #leaving this as an arbitrary lower bound on K

    #if params dict has a K_cap parameter, use this to set the max cell value for the K raster
    if 'K_cap' in params.keys():
        K[K>float(params['K_cap'])] = float(params['K_cap'])
        print(K.max())
    assert K.min() >= 0
    assert True not in np.isnan(K)
    assert True not in np.isinf(K)


    #print('Plotting K')
    #plt_K = landscape.Landscape(land.dims, K)
    #pop.show(plt_K)
    #raw_input()


    #print('Plotting n_pairs')
    #plt_n_pairs = landscape.Landscape(land.dims, n_pairs)
    #pop.show(plt_n_pairs)
    #raw_input()



    #use N, K, and params['r'] to calc dN/dt raster, according to the chosen pop growth eqxn
        #NOTE: Currently, only option and default is basic logistic eqxn: dN/dt = r*(1-N/K)*N
    dNdt = calc_dNdt(land, N, K, params, pop_growth_eq = 'logistic').raster
    dNdt[N == 0] = 0
    assert True not in np.isnan(dNdt)
    assert True not in np.isinf(dNdt), 'The following cells are infinite: \n\t%s' % str([i for i, n in enumerate(dNdt.ravel()) if np.isinf(n)])


    #print('Plotting dNdt')
    #plt_dNdt = landscape.Landscape(land.dims, dNdt)
    #pop.show(plt_dNdt)
    #raw_input()



    #use n_pairs, params['b'] (i.e. probability that a pair gives birth), and params['lambda_offspring'] (i.e.
    #expected number of births per pair) to calculate the N_b raster (expected number of births in each cell)
    #according to N_b = b*E[N_b/pair]*n_pairs where E[N_b/pair] = lambda (b/c births Poission distributed)
        #NOTE: Would be good to in the future to use neg_binom dist to allow births to be modeled as Poisson
        #(using correct parameterization of neg_binom to approximate that) or overdispersed (to be parmeterized with field data)

        #NOTE: Does it matter that I'm not attempting to account for the fact that pairs will disperse out of
        #the cell? Almost certainly a larger effect with increasing dispersal-kernel distance. But nonetheless,
        #would it even make sense (and even BE TRACTABLE) to calculate it any other way? My intuition is NO, but worth considering
        #further later on...

        #NOTE: Currently, birth rate is intrinsic to the species and cannot be made a function of density or
        #spatial fitness. Might be worth adding this at some point in the future though.
    N_b = params['b'] * params['lambda_offspring'] * n_pairs
    assert N_b.min() >= 0
    assert True not in np.isnan(N_b)
    assert True not in np.isinf(N_b)

   
    #print('Plotting N_b')
    #plt_N_b = landscape.Landscape(land.dims, N_b)
    #pop.show(plt_N_b)
    #raw_input()



    #Use N_b - dN/dt to calculate N_d raster (expected number of deaths in each cell)
    N_d = N_b - dNdt
    assert True not in np.isnan(N_d)
    assert True not in np.isinf(N_d)

   
    #print('Plotting N_d')
    #plt_N_d = landscape.Landscape(land.dims, N_d)
    #pop.show(plt_N_d)
    #raw_input()



    #Use N_d/N to calculate d raster (probability of death of individuals in each cell)
    #d = N_d/N
    #d[d<0] = 0

    #Instead of the above, just normalize N_d to use it as the raster of the probability of death of individuals in each cell
    d = (N_d - N_d.min())/(N_d.max() - N_d.min())
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
    d[d<d_min] = d_min
    d[d>d_max] = d_max


    assert d.min() >= 0, 'd.min() is %0.2f' % d.min()
    assert d.max() <= 1, 'd.max() is %0.2f' % d.max()
    assert True not in np.isnan(d)
    assert True not in np.isinf(d)


    #print('Plotting d')
    #print('Max val: %0.2f' % d.max())
    #plt_d = landscape.Landscape(land.dims, d)
    #pop.show(plt_d)
    #raw_input()


    #Now implement births
        #Feed the mating pairs and params['b'] to the mating functions, to produce and disperse zygotes
    pop.mate(land, params, pairs)

    #Now implement deaths
        #NOTE: The way I'm looking at it now, calculating likely births from number of existing pairs, then
        #calculating death raster, then having births take place, then having per-individual death probabilies
        #determined and having deaths take place (such that newborns could also 
        #immediately die after birth, based on the death raster calculated before they were born, and even
        #based on individs who weren't their parents, assuming they dispersed across cell lines), which is what
        #I'm doing currently, is the only obvious/reasonable way to do this. Of course there will always be an
        #inherent rub between reality and discrete-time simulations because of the continuous/discrete time
        #conflict. But in this case, is there any strong reason to CHANGE THIS? Or are there any foreseeable
        #and undesirable results/effects?... NEED TO PUT MORE THOUGHT INTO THIS LATER.

    death_probs = dict([(i, d[int(ind.y), int(ind.x)]) for i,ind in pop.individs.items()])

    #plt.hist(death_probs.values(), bins = 50)
    #raw_input()

    #If selection = True, then use the d raster and individuals' relative fitnesses to calculate
    #per-individual probabilities of death
    if selection == True:
        pass
    
    #If age_stage_d <> None then use the age_stage_d array, and pop.get_age(), and the d raster 
    #to calculate per-individual death probabilities
    if age_stage_d <> None:
        pass
        
    #Feed the per-individual death probabilities into the deaths function, to generate deaths and cull
    #those individuals
    num_deaths = mort(land, pop, params, death_probs)

    print '\n\t%i individuals dead' % num_deaths

    #Check if extinct
        #If so, will return that information and exit
    pop.check_extinct()

    #Age the population NOTE: SHOULD THIS HAPPEN LATER, SO THAT NEWBORN ACTUALLY LIVE THROUGH ONE FULL CYCLE?
    #pop.birthday()

    #If burn_in == True:
    if burn_in == True:
        pass
        #Update pop.K
        #Run the burn_in functions to update the running metrics, and assess stationarity, and return decision
        #of whether or not stationarity has been achieved

        #return decision
        
    #Else exit
    else:
        pass




###############################
#NOTE: MOVE TO burn_in.py

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



