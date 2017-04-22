#!/usr/bin/python
#PRACTICE_burn_in.py

'''Defines the burn-in function.'''


def burn_in(pop, land, params): 
    print '\n\nSTARTING BURN-IN.\n\t(Will run for %i timesteps.)\n\n' % params['burn_T']
    #pop.show(land = land, colorbar = True)

    #add functionality for running multiple pops at once
    multi_pop = False
    if type(pop) in (list, tuple):
        multi_pop = True


    if multi_pop == False:
        pop.set_K(pop.calc_density(land = land, set_N = False))
        for burn_t in range(params['burn_T']):
            print('Timestep %i:' % burn_t)
            pop.birthday()
            pop.set_K(pop.calc_density(land = land, set_N = False, max_val = params['K_cap']))
            pop.mutate(params = params, t = burn_t)
            pop.move(land = land, params = params)
            pop.calc_density(land = land, window_width = max(1, params['mu_distance']), set_N = True)
            demography.pop_dynamics(land = land, pop = pop, params = params, selection = False, burn_in = True, d_min = 0.1, d_max = 0.9)
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


