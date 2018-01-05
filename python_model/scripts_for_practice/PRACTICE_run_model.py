execfile('./params.py')
execfile('./scripts_for_practice/PRACTICE_create_land_genomic_arch_pop.py')
execfile('./scripts_for_practice/PRACTICE_imports_and_reloads.py')

pop.set_K(land.scapes[params['movement_surf_scape_num']])

break_burn_in = False
burn_in_test_t = params['burn_T_min']
t = 0
while break_burn_in == False:
	print('###############\n\n TIMESTEP %i\n' % t)
	pop.increment_age_stage()
	pop.mutate(params,t);pop.set_Nt()
	pop.move(land, params) 
	extinct = demography.pop_dynamics(land, pop, params, with_selection = False) 
	t +=1
	if extinct == 1: 
		break 
	break_burn_in = (len(pop.Nt) > burn_in_test_t and burn_in.adf_threshold_test(pop, burn_in_test_t, 0.05) and burn_in.tt_threshold_test(pop, burn_in_test_t, 0.05))


print('~~~~~~~~~~~~~~~~\n\n\n\t\tBURN-IN COMPLETE\n\n\n~~~~~~~~~~~~~~~~~~')
    

for t in range(params['T']):
	print('###############\n\n TIMESTEP %i\n' % t)
	pop.increment_age_stage()
	pop.mutate(params,t);pop.set_Nt()
	pop.move(land, params)
	extinct = demography.pop_dynamics(land, pop, params, with_selection = True) 
	if extinct == 1: 
		break

