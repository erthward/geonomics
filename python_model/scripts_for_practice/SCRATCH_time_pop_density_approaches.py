import time
import numpy.random as r

hwws = [1, 2.5, 5, 10, 20]

ltys = ['-', '--', ':']

fig = plt.figure()

fig.suptitle('Comparing 2 pop-density calculation approaches (red and blue)\nby landscape dims, pop size (500,1500,5000: solid, dashed, dotted), and neighborhood count window widths')

for n, d in enumerate([50,100,200]):
    land = landscape.Landscape((d,d), zeros((d,d)))
    ax = fig.add_subplot(2,2,n+1)
    ax.set_title('LANDSCAPE: %ix%i' % (d, d))
    ax.set_xlabel('half-window width (cells)')
    ax.set_ylabel('time (sec)')

    for m, p in enumerate([500,1500,5000]):
            xs = r.uniform(0,d,p)
            ys = r.uniform(0,d,p)
            pop = population.Population(p, dict([(i, individual.Individual(np.array([0]), xs[i], ys[i])) for i in range(len(xs))]), g, p, 100)
            alt_t = []
            alt2_t = []
            for hww in hwws:
                search_area_array = create_search_area_array(hww, (d,d))
                alt_times = []
                alt2_times = []
                for i in range(5):
                    start = time.time() 
                    alt_calc_density(pop, land)
                    stop = time.time()
                    alt_times.append(stop-start)

                    start = time.time()
                    alt2_calc_density(pop, land, search_area_array, window_width = 2*hww)
                    stop = time.time()
                    alt2_times.append(stop-start)

                alt_t.append(mean(alt_times))
                alt2_t.append(mean(alt2_times))
        
            plt.plot(hwws, alt_t, 'r%s' % ltys[m])
            plt.plot(hwws, alt2_t, 'b%s' % ltys[m])

plt.show()
