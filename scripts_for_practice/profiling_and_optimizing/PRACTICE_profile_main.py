def __profile_main__():
    pop.increment_age_stage(burn = False)
    pop.set_Nt()
    pop.move(land, params)
    extinct = demography.pop_dynamics(land, pop, params, with_selection = True) 
    pop.mutate(log = False)
    return

