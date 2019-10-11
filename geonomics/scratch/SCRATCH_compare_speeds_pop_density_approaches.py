import time
import numpy.random as r

wws = [1, 2.5, 5, 10, 20]

ltys = ['-', '--', ':']

fig = plt.figure()

fig.suptitle('Comparing 2 pop-density calculation approaches (original:red, string-based cell matching:blue)\nby landscape dims, pop size (500,1500,5000: solid, dashed, dotted), and neighborhood count window widths')

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
            new_t = []
            for q, ww in enumerate(wws):
                search_area_array = create_search_area_array(ww, (d,d))
                alt_times = []
                new_times = []
                for i in range(5):
                    start = time.time() 
                    pop.calc_density(land, window_width = ww)
                    stop = time.time()
                    alt_times.append(stop-start)

                    start = time.time()
                    gs = Grid_Stack(land, ww)
                    gs.calc_density(list(pop.get_x_coords().values()), list(pop.get_y_coords().values()))
                    stop = time.time()
                    new_times.append(stop-start)

                alt_t.append(mean(alt_times))
                new_t.append(mean(new_times))
        
            plt.plot(wws, alt_t, 'r%s' % ltys[m])
            plt.plot(wws, new_t, 'b%s' % ltys[m])

fig.savefig('Users/mariko_terasaki/Desktop/drew/geonomics_presentations/compare_speeds_pop_density_algorithms.pdf')
plt.show()
