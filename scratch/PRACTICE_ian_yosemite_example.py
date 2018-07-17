exec(open('./params.py', 'r').read())
exec(open('./scratch/PRACTICE_imports_and_reloads.py', 'r').read())
exec(open('./scratch/PRACTICE_create_land_genomic_arch_pop.py', 'r').read())

run = True

#little fn to calculate a habitat raster that's 1 at the center of env-var's range and drops toward 0 at extremes
def calc_hab_rast(env_rast):
    hab = 1-abs(env_rast-0.5)
    hab = spt.scale_raster(hab)[0]
    return(hab)

#read in and reload scape #1 (because I must be already 0-1-scaling it in one of the other scripts that have been called by this point... NOTE: investigate this)
rast, dim, res, ulc = io.read_raster('/home/ihavehands/Desktop/yos_30yr_normals_90x90.tif')
land.scapes[1].raster = rast

#extract the historical raster (read in from yos_30yr_normals_90x90.tif)
hist_rast = np.flipud(land.scapes[1].raster)

#and calculate a future-change raster, where change is fastest at highest elevations
#(add 2 Kelvin raster wide plus an additional fraction of 2 that is greatest at highest elevations)
fut_rast = hist_rast + 2 + 2 * (hist_rast.max()-hist_rast)/(hist_rast.max()-hist_rast.min())

#now scale both to 0<=x<=1, using the min and max across the changing time-series
sc_hist_rast, hist_min, hist_max = spt.scale_raster(hist_rast, hist_rast.min(), fut_rast.max())
sc_fut_rast, fut_min, fut_max = spt.scale_raster(fut_rast, hist_rast.min(), fut_rast.max())

#set raster 1 to the scaled historical raster
land.scapes[1].raster = sc_hist_rast

#and set scape 0 to the habitat raster calculated from this
land.scapes[0].raster = calc_hab_rast(sc_hist_rast)

#and set the time-series end-scapes 
params['change']['land'][1]['end_scape'] = sc_fut_rast
params['change']['land'][0]['end_scape'] = calc_hab_rast(sc_fut_rast)

#recalculate the starting movement surface
land.movement_surf = spt.Movement_Surface(land)

#set the carrying capacity to 3, then a weighted version of that raster
pop.set_K(land.scapes[0].raster*3)
pop.K = pop.K*0.005 + pop.K*(0.995*(pop.K-pop.K.min()))*0.5
#pop.K = pop.K*5

#then create the Land_Changer
lc = change.Land_Changer(land, params)

#tweak the genomic architecture
ga = pop.genomic_arch
t = ga.traits[0]
ga.p[t.loci] = 0.5
t.alpha = np.array([0.1,-0.1]*5)

def burn(stop_after = None):
    break_burn_in = False
    burn_in_test_t = params['model']['burn']['burn_T_min']
    t = 0
    while break_burn_in == False:
        print('###############\n\n TIMESTEP %i' % t)
        print('     POP %i\n' % pop.get_size())
        pop.reset_age_stage(burn=True)
        pop.set_Nt()
        pop.do_movement(land)
        extinct = demography.do_pop_dynamics(land, pop, with_selection=False, burn=True)
        t += 1
        if extinct == 1:
            break
        break_burn_in = (len(pop.Nt) > burn_in_test_t and burn_in.test_adf_threshold(pop, burn_in_test_t, 0.05) and burn_in.test_tt_threshold( pop, burn_in_test_t, 0.05))
        if stop_after is not None and t == stop_after:
            break


    print('~~~~~~~~~~~~~~~~\n\n\n\t\tBURN-IN COMPLETE\n\n\n~~~~~~~~~~~~~~~~~~')


def main(T, reassign_genomes=False, fits = None):
    if reassign_genomes == True:
        print('\n\nReassigning genomes...\n\n')
        genome.reset_genomes(pop, params)
        [i.set_phenotype(genomic_arch) for i in pop.individs.values()];
    for t in range(T):
        if t == 1399:
            ax = fig.add_subplot(232)
            ax.set_title('Timestep 1400')
            pop.show_phenotype(0, land, scape_num = 0, size = 7)
            ax = fig.add_subplot(235)
            pop.show_phenotype(0, land, scape_num = 1, size = 7)
        if t == 1999:
            ax = fig.add_subplot(233)
            ax.set_title('Timestep 2000')
            pop.show_phenotype(0, land, scape_num = 0, size = 7)
            ax = fig.add_subplot(236)
            pop.show_phenotype(0, land, scape_num = 1, size = 7)
        print('###############\n\n TIMESTEP %i' % t)
        print('     POP %i\n' % pop.get_size())
        lc.make_change(pop.t)
        pop.reset_age_stage(burn=False)
        pop.set_Nt()
        pop.do_movement(land)
        extinct = demography.do_pop_dynamics(land, pop, with_selection=True)
        if extinct == 1:
            break
        fits.append(mean(pop.get_fitness()))
        pop.do_mutation(log=False)

if run:
    fig = plt.figure()
    plt.suptitle('Yosemite time series:\naltitude-weighted 2-to-4 K increase in 30-yr. normals\nselection by temperature, 10-genic trait')
    ax = fig.add_subplot(231)
    ax.set_title('Timestep 0')
    burn()
    main(0, True)
    pop.show_phenotype(0, land, scape_num = 0, size = 7)
    ax = fig.add_subplot(234)
    pop.show_phenotype(0, land, scape_num = 1, size = 7)
    fits = []
    main(2000, False, fits)
    plt.show()
   
    fig2 = plt.figure()
    plt.suptitle('Yosemite example, mean fitness as a function of time')
    plt.plot(range(len(fits)), fits)
    plt.xlabel('timesteps')
    plt.ylabel('fitness')
    plt.show()

