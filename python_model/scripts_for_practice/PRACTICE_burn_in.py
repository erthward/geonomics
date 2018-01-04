#!/usr/bin/python
#PRACTICE_burn_in.py

'''Defines the burn-in function.'''

from collections import Counter as C

def burn_in(pop, land, params, het, maf): 
    print '\n\nSTARTING BURN-IN.\n\t(Will run for %i timesteps.)\n\n' % params['burn_T']
    #pop.show(land = land, colorbar = True)

    #add functionality for running multiple pops at once
    multi_pop = False
    if type(pop) in (list, tuple):
        multi_pop = True


    if multi_pop == False:
        #set starting K_loc to the starting pop density
        #pop.set_K(pop.calc_density(land = land, set_N = False))
        #NOTE: Instead, just setting it statically to K_cap * the movement raster (which can be seen simultaneously as a 'hab qual' raster)
        import landscape 
        K = landscape.Landscape(land.dims, params['K_cap']*land.scapes[params['movement_surf_scape_num']].raster)
        #K.raster[K.raster<0.5] = 0.5
        pop.set_K(K)
        for burn_t in range(params['burn_T']):
            print('Timestep %i:' % burn_t)
            cts = C(pop.get_genotype(0,0).values())
            het.append(cts[0.5]/float(pop.census()))
            maf.append((cts[1]*2 + cts[0.5])/(2.*pop.census()))
            #print('\tbirthday\n')
            pop.birthday()
            #print('\tset_K\n')
            #NOTE: updating K was leading to runaway pop growth, so arresting it for the moment
            #pop.set_K(pop.calc_density(land = land, set_N = False, max_val = params['K_cap']))
            #print('\tmutate\n')
            pop.mutate(params = params, t = burn_t)
            #print('\tmove\n')
            pop.move(land = land, params = params)
            #print('\tcalc_density\n')
            #pop.calc_density(land = land, set_N = True)#window_width = max(1.01, params['mu_distance']), set_N = True)
            #print('\tpop_dynamics\n')
            if burn_t < 500:
                demography.pop_dynamics(land = land, pop = pop, params = params, selection = False, burn_in = True)
                #NOTE: CAN SWITCH TO THIS LINE TO START POP-DYN DEBUGGING BEHAVIOR AT CERTAIN TIMESTEP; FEED TIMESTEP IN FOR PLOT TITLES
                #demography.pop_dynamics(land = land, pop = pop, params = params, selection = False, burn_in = True, debug = burn_t)
            else:
                demography.pop_dynamics(land = land, pop = pop, params = params, selection = True, burn_in = True)

            #NOTE: Making the d_min and d_max values considerably more permissive than the default settings, just for the burn-in period, to allow for more pronounced shifts in spatial distribution of individs during the iterative algorithm without too much 'penalty'

            #pop.select(t = burn_t, params = params)
            #pop.check_extinct()
            print '\n\n%i timesteps run.  Current status:\n\n\t%i individuals\n\n' % (burn_t+1, pop.census())
            print '\n----------------------------------------------------------\n\n' 

        print '\n\nBURN-IN FINISHED.\n\n'
        #pop.show(land = land, colorbar = False)
        #mpl.pyplot.close()


    #NOTE: NOT KEPT UP TO DATE WITH THE VERSION ABOVE!
    elif multi_pop == True:
        for burn_t in range(params['burn_T']):
            print('Timestep %i:' % burn_t+1)
            n = 1
            for p in pop:
                print('\tPopulation %i:' % n)
                p.birthday()
                p.move(land = land, params = params)
                p.mate(land = land, params = params)
                p.select(t = burn_t, params = params)
                p.mutate(params = params, t = burn_t)
                p.check_extinct()
                n += 1
            print('\n\n%i timesteps run.   Current statuses:\n\t%s' % ( burn_t + 1, '\n\t'.join(['\t'.join(['POP%i' % (n+1), '%i individuals' % p.census()]) for n, p in enumerate(pop)]))) 


            print '\n----------------------------------------------------------\n\n' 


        print '\n\nBURN-IN FINISHED.\n\n'
        #pop.show(land = land, colorbar = False)
        #mpl.pyplot.close()

