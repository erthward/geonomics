#!/usr/bin/python

from time import time
import matplotlib.pyplot as plt

w = [1,2,3,4,5,6,7,8,9,10,15,20,25,35,50]

times = []

for val in w:
    tots = []
    for it in range(20):
        start = time()
        pop.calc_density(land, window_width = val)
        stop = time()
        tots.append(stop-start)
    times.append(mean(tots))


plt.plot(w, times)


