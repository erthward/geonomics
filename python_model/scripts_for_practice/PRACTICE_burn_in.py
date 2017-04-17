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
        for burn_t in range(params['burn_T']):

            print('Timestep %i:' %burn_t)
            pop.birthday()
            pop.move(land = land, params = params)
            pop.mate(land = land, params = params)
            pop.select(t = burn_t, params = params)
            pop.mutate(params = params, t = burn_t)
            pop.check_extinct()
            print '\n\n%i timesteps run.  Current status:\n\n\t%i individuals\n\n' % (burn_t+1, pop.census())
            print '\n----------------------------------------------------------\n\n' 

        print '\n\nBURN-IN FINISHED.\n\n'
        #pop.show(land = land, colorbar = False)
        #mpl.pyplot.close()



    elif multi_pop == True:
        for burn_t in range(params['burn_T']):
            print('Timestep %i:' %burn_t)
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


