import numpy.random as rand

def reassign_coords(pop):
    while 0.0 in [v[1] for v in pop.get_habitat().values()]:
        inds = [i for i,v in pop.get_habitat().items() if v[1] == 0.0]
        print(inds)

        new_x = rand.uniform(low = 0, high = 50, size = len(inds))
        new_y = rand.uniform(low = 0, high = 50, size = len(inds))
        print(new_x, new_y)

        for i in range(len(inds)):
            pop.individs[inds[i]].x = new_x[i] 
            pop.individs[inds[i]].y = new_y[i]
