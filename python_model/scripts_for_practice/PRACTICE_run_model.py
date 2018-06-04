exec(open('./params.py', 'r').read())
exec(open('./scripts_for_practice/PRACTICE_imports_and_reloads.py', 'r').read())
exec(open('./scripts_for_practice/PRACTICE_create_land_genomic_arch_pop.py', 'r').read())

pop.set_K(land.scapes[params['n_movement_surf_scape']].raster)

pop.genomic_arch.traits[0].s = 0.1

# pop.genomic_arch.traits[0].alpha = np.array([0.25])

for i in pop.individs.keys():
    pop.individs[i].genome[pop.genomic_arch.traits[0].loci, :] = r.binomial(1, 0.5, 2)
    

def burn(stop_after = None):
    break_burn_in = False
    burn_in_test_t = params['burn_T_min']
    t = 0
    while break_burn_in == False:
        print('###############\n\n TIMESTEP %i' % t)
        print('     POP %i\n' % pop.census())
        pop.increment_age_stage(burn=True)
        # pop.mutate(params,t)
        pop.set_Nt()
        pop.move(land, params)
        extinct = demography.pop_dynamics(land, pop, params, with_selection=False, burn=True)
        t += 1
        if extinct == 1:
            break
        break_burn_in = (len(pop.Nt) > burn_in_test_t and burn_in.adf_threshold_test(pop, burn_in_test_t, 0.05) and burn_in.tt_threshold_test( pop, burn_in_test_t, 0.05))
        if stop_after is not None and t == stop_after:
            break


    print('~~~~~~~~~~~~~~~~\n\n\n\t\tBURN-IN COMPLETE\n\n\n~~~~~~~~~~~~~~~~~~')


def main(T, reassign_genomes=False):
    if reassign_genomes == True:
        print('\n\nReassigning genomes...\n\n')
        genome.reassign_genomes(pop, params)
        [i.set_phenotype(genomic_arch) for i in pop.individs.values()];
    for t in range(T):
        print('###############\n\n TIMESTEP %i' % t)
        print('     POP %i\n' % pop.census())
        pop.increment_age_stage(burn=False)
        pop.set_Nt()
        pop.move(land, params)
        extinct = demography.pop_dynamics(land, pop, params, with_selection=True)
        if extinct == 1:
            break
        pop.mutate(log=False)

