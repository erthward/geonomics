import geonomics as gnx

from copy import deepcopy

orig_params = gnx.read_parameters_file('./test/validation/wf/wf_params_DEH.py')

K_factors = [10, 20, 50]
mean_t_fix_list = []
mean_pop_size_list = []

def get_allele_freqs(spp):
    freqs = {}
    for loc in range(spp.gen_arch.L):
        freq = sum(np.hstack([ind.genome[loc,:] for ind in spp.values(
            )]))/(2*len(spp))
        freqs[loc] = freq
    return freqs

for K_fact in K_factors:
    params = deepcopy(orig_params)
    params.comm.species['spp_0'].init['K_factor'] = K_fact
    print("USING K_fact %0.2f" % params.comm.species['spp_0'].init.K_factor)

    mod = gnx.make_model(params)

    mod.walk(mode = 'burn', T = 10000, verbose = True)

    freqs = {loc:[] for loc in range(mod.comm[0].gen_arch.L)}
    #run model until all alleles have fixed
    while (False not in [len(f) == 0 for f in freqs.values()]
        or False in [f[-1] in (0,1) for f in freqs.values()]):
    #for t in range(mod.params.model.T):
        freqs_t = get_allele_freqs(mod.comm[0])
        for loc, freq in freqs_t.items():
            freqs[loc].append(freq)
        mod.walk(mode = 'main', T = 1, verbose = True)
        print('\n%i ALLELES FIXED\n' % sum(
            [f[-1] in (0, 1) for f in freqs.values()]))


    #calculate average time to fixation and average pop size, and record
    mean_pop_size = mean(mod.comm[0].Nt)
    t_fix_list = [f.index(0) if 0 in f else f.index(1) for f in freqs.values()]
    mean_t_fix = mean(t_fix_list)
    mean_pop_size_list.append(mean_pop_size)
    mean_t_fix_list.append(mean_t_fix)

    fig = plt.figure()
    plt.suptitle(("Drift trajectories for %i independent loci\nin Wright-Fisher "
        "approximation validation test, mean pop. size = %0.2f") % (
        mod.comm[0].gen_arch.L, mean_pop_size))
    plt.xlabel('model time (timesteps)')
    plt.ylabel('frequency of 1-allele')
    plt.xlim(0,len(mod.comm[0].Nt))
    plt.ylim(0,1)
    for loc, freq_list in freqs.items():
        plt.plot(range(len(freq_list)), freq_list, '-')
    plt.show()

fig = plt.figure()
plt.suptitle(("Mean population size and mean time to fixation,\nas a function "
    "of K_factor, for %i independent, neutral loci") % mod.comm[0].gen_arch.L)
ax = fig.add_subplot(1,2,1)
plt.title('Mean population size')
plt.xlabel('K_factor')
plt.ylabel('mean population size')
plt.plot(K_factors, mean_pop_size_list)
ax = fig.add_subplot(1,2,2)
plt.title('Mean time to fixation')
plt.xlabel('K_factor')
plt.ylabel('mean time to fixation')
#plot the data
plt.plot(K_factors, mean_t_fix_list)
#plot expected mean fixation times according to W-F model
#(according to theory, expected time to fixation for loci starting at 0.5/0.5
#allele freqs is 2.776N, i.e. -4N(log(p)p + log(q)q) solved for p = q = 0.5;
#here I expect some deviation from that because this
#is based on a single simulated instantiation, and also because my
#populations are approximating W-F populations but not perfectly, so most
#likely N_t > N_e)
plt.plot(K_factors, [2.776*n for n in mean_pop_size_list], 'or')

