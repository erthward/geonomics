#exec(open('./params.py', 'r').read())
#exec(open('./scratch/PRACTICE_imports_and_reloads.py', 'r').read())
#exec(open('./scratch/PRACTICE_create_land_genomic_arch_pop.py', 'r').read())

#pop.set_K(land[params['land']['movement_surf']['movement_surf_scape_num']].raster)

#lc = Land_Changer(land, params)

def burn(stop_after = None, burn = True):
    break_burn_in = False
    burn_in_test_t = params.model.time.burn_T
    t = 0
    while break_burn_in == False:
        print('###############\n\n TIMESTEP %i' % t)
        print('     POP %i\n' % len(pop))
        pop.set_age_stage(burn=True)
        pop.set_Nt()
        pop.do_movement()
        extinct = demography.do_pop_dynamics(land, pop, with_selection=False, burn=True)
        t += 1
        if extinct == 1:
            break
        break_burn_in = (len(pop.Nt) > burn_in_test_t and burnin.test_adf_threshold(pop, burn_in_test_t, 0.05) and burnin.test_t_threshold( pop, burn_in_test_t, 0.05))
        if stop_after is not None and t == stop_after:
            break


    print('~~~~~~~~~~~~~~~~\n\n\n\t\tBURN-IN COMPLETE\n\n\n~~~~~~~~~~~~~~~~~~')


def main(T, reassign_genomes=True):
    if reassign_genomes == True:
        print('\n\nReassigning genomes...\n\n')
        genome.set_genomes(pop, params.model.time.burn_T, 1000)
        [ind.set_phenotype(pop.gen_arch) for ind in pop.inds];
    for t in range(T):
        print('###############\n\n TIMESTEP %i' % t)
        print('     POP %i\n' % len(pop))
        #lc.make_change(t)
        pop.set_age_stage(burn=False)
        pop.set_Nt()
        pop.do_movement()
        extinct = demography.do_pop_dynamics(land, pop, with_selection=True)
        if extinct == 1:
            break
