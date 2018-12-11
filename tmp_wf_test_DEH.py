import geonomics as gnx

mod = gnx.make_model('./test/validation/wf/wf_params_DEH.py')

mod.walk(mode = 'burn', T = 10000, verbose = True)

def get_allele_freqs(spp):
    freqs = {}
    for loc in range(spp.gen_arch.L):
        freq = sum(np.hstack([ind.genome[loc,:] for ind in spp.values(
            )]))/(2*len(spp))
        freqs[loc] = freq
    return freqs

freqs = {loc:[] for loc in range(mod.comm[0].gen_arch.L)}
for t in range(mod.params.model.T):
    freqs_t = get_allele_freqs(mod.comm[0])
    for loc, freq in freqs_t.items():
        freqs[loc].append(freq)
    mod.walk(mode = 'main', T = 1, verbose = True)


fig = plt.figure()
plt.suptitle(("Drift trajectories for 100 independent loci\nin Wright-Fisher "
    "approximation validation test."))
plt.xlabel('model time (timesteps)')
plt.ylabel('frequency of 1-allele')
plt.xlim(0,500)
plt.ylim(0,1)
for loc, freq_list in freqs.items():
    plt.plot(range(len(freq_list)), freq_list, '-')
plt.show()

