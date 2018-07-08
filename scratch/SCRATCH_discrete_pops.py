#!/usr/bin/python

#assuming PRACTICE_run_model.py has already been run, but burn-in hasn't yet happened:

n_pops = 5
pop_buf = 2
K_factor = 20




new_land = np.zeros(land.dims)

pop_locs = r.choice(range(land.dims[0]), 2*n_pops).reshape((n_pops, 2))

for pop_loc in pop_locs:
    new_land[pop_loc[0]-pop_buf:pop_loc[0]+pop_buf, pop_loc[1]-pop_buf:pop_loc[1]+pop_buf] = 1

land.scapes[params['n_movement_surf_scape']].raster = new_land

pop.K = K_factor*new_land

land.movement_surf = movement.create_movement_surface(land, params)


