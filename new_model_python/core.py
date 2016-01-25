#!/usr/bin/python
#core.py

                                                                                                                                                         
#---------#
# BURN IN #
#---------#


def burn_in(land, pop, params):

    print '\n\nSTARTING BURN-IN.\n\t(Will run for %i timesteps.)\n\n' % params['burn_T']

    pop.show(land = land, colorbar = True)

    for burn_t in range(params['burn_T']):
    
        pop.birthday()
    
        pop.move(land = land, params = params)
    
        pop.mate(land = land, params = params)

        pop.select(t = burn_t, params = params)
    
        pop.mutate(params = params, t = burn_t)
    
        pop.check_extinct()

        
        print '\n\n%i timesteps run. Current status:\n\n\t%i individuals\n\n' % (burn_t+1, pop.census())

        print '\n----------------------------------------------------------\n\n'

    
    print '\n\nBURN-IN FINISHED.\n\n'

    pop.show(land = land, colorbar = False)

    #mpl.pyplot.close()

