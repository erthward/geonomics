def old_main(T, reassign_genomes=False):
    if reassign_genomes == True:
        print('\n\nReassigning genomes...\n\n')
        genome.reassign_genomes(pop, params)
    for t in range(T):
        print('###############\n\n TIMESTEP %i\n' % t)
        pop.increment_age_stage(burn=False)
        pop.set_Nt()
        pop.move(land, params)
        extinct = demography.pop_dynamics(land, pop, params, with_selection=True)
        if extinct == 1:
            break
        pop.mutate(log=False)



def alt_main(T, reassign_genomes=False):
    if reassign_genomes == True:
        print('\n\nReassigning genomes...\n\n')
        genome.reassign_genomes(pop, params)
    for t in range(T):
        print('###############\n\n TIMESTEP %i\n' % t)
        pop.increment_age_stage(burn=False)
        pop.set_Nt()
        #pop.move(land, params)
        alt_move(pop, land, params)
        pop.set_habitat(land)
        extinct = demography.pop_dynamics(land, pop, params, with_selection=True)
        if extinct == 1:
            break
        pop.mutate(log=False)



def alt_alt_main(T, reassign_genomes=False):
    if reassign_genomes == True:
        print('\n\nReassigning genomes...\n\n')
        genome.reassign_genomes(pop, params)
    for t in range(T):
        print('###############\n\n TIMESTEP %i\n' % t)
        pop.increment_age_stage(burn=False)
        pop.set_Nt()
        #pop.move(land, params)
        alt_alt_move(pop, land, params)
        pop.set_habitat(land)
        extinct = demography.pop_dynamics(land, pop, params, with_selection=True)
        if extinct == 1:
            break
        pop.mutate(log=False)





