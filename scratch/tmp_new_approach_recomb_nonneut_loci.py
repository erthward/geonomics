import numpy as np
import matplotlib.pyplot as plt
import bisect
import time

n_loci_list = [10, 100]#, 1000, 10000]
n_gametes_list = [10, 100]#, 1000, 10000]


def time_recombs(n_loci, n_gametes):
    # non-neutral loci
    locs = list(set(np.random.randint(low=0, high=10000-1, size = 10500)))[:n_loci] 
    locs.sort()

    # recombiation event breakpoints (minus the 0 and max-length positions)
    recombs = []
    for i in range(n_gametes):
        n_recombs = np.random.poisson(100)
        recomb = np.random.choice(range(n_loci), n_recombs, replace=False)
        recomb.sort()
        recombs.append(recomb)

    # non-neutral genotypes for 10_000 individuals
    gts = [np.random.binomial(1, 0.5, size=2*len(locs)).reshape(
                    (2, len(locs))) for _ in range(len(recombs))]

    start = time.time()
    # now get the gametes produced from each of those genotypes
    res = [gts[n][[(np.random.binomial(
                    0, 0.5) + (bisect.bisect_right(
                    r, l) % 2))%2 for l in locs], range(len(
                    locs))] for n, r in enumerate(recombs)]
    stop = time.time()
    diff = stop - start
    print('This took %0.2f seconds.' % (diff))
    return(diff)

# time them
loci_times = {n_loci: [] for n_loci in n_loci_list}
gam_times = {n_gam: [] for n_gam in n_gametes_list}
for n in n_loci_list:
    for m in n_gametes_list:
        print('n loci', n)
        print('n gams', m)
        t = time_recombs(n, m)
        loci_times[n].append(t)
        gam_times[m].append(t)

print(loci_times)
print(gam_times)

# plot them
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax1.set_title("loci times")
ax1.set_xlabel("n loci")
ax1.set_ylabel("time")
ax2 = fig.add_subplot(122)
ax2.set_title("loci times")
ax2.set_xlabel("n gams")
ax2.set_ylabel("time")

for i, item in enumerate(loci_times.items()):
    print(n_gametes_list)
    print(item[1])
    ax2.plot(n_gametes_list, item[1], label=str(item[0]))
lg2 = ax2.legend()
lg2.set_title('n loci')

for i, item in enumerate(gam_times.items()):
    ax1.plot(n_loci_list, item[1], label=str(item[0]))
lg1 = ax1.legend()
lg1.set_title('n gams')
plt.show()


# NEXT: then would gather those gametes in twos to produce new offspring genomes,
# and use the genomes to calculate phenotypes and fitnesses for each offspring
