import time
import numpy.random as r

hwws = [1, 2.5, 5, 10, 20]
bws = [0.02*n for n in hwws]

ltys = ['-', '--', ':']

fig = plt.figure()

fig.suptitle('Comparing 2 pop-density calculation approaches (original:red, KDE:blue)\nby landscape dims, pop size (500,1500,5000: solid, dashed, dotted), and neighborhood count window widths')

for n, d in enumerate([50,100,200]):
    land = landscape.Landscape((d,d), zeros((d,d)))
    ax = fig.add_subplot(2,2,n+1)
    ax.set_title('LANDSCAPE: %ix%i' % (d, d))
    ax.set_xlabel('half-window width (cells)')
    ax.set_ylabel('time (sec)')

    for m, p in enumerate([500,1500,5000]):
            xs = r.uniform(0,d,p)
            ys = r.uniform(0,d,p)
            pop = population.Population(p, dict([(i, individual.Individual(np.array([0]), xs[i], ys[i])) for i in range(len(xs))]), ga, p, 100)
            pop.Nt = [p]
            alt_t = []
            kde_t = []
            for q, hww in enumerate(hwws):
                search_area_array = create_search_area_array(hww, (d,d))
                alt_times = []
                kde_times = []
                for i in range(5):
                    start = time.time() 
                    alt_calc_density(pop, land)
                    stop = time.time()
                    alt_times.append(stop-start)

                    start = time.time()
                    calc_kde_pop_density(pop, land, bws[q])
                    stop = time.time()
                    kde_times.append(stop-start)

                alt_t.append(mean(alt_times))
                kde_t.append(mean(kde_times))
        
            plt.plot(hwws, alt_t, 'r%s' % ltys[m])
            plt.plot(hwws, kde_t, 'b%s' % ltys[m])

plt.show()
