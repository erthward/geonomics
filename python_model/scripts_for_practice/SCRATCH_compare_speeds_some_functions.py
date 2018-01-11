#!/usr/bin/python


import numpy as np
import numpy.random as r
import matplotlib as mpl
import matplotlib.pyplot as plt
import timeit


#NOTE: I will probably want to look into Cython for some of the most time-consuming operations!


###############
# RECOMBINATION
###############


n_loci = [100,1000,10000]

#create recomb_vectors of diff lengths:
    #one to use in the current approach
curr_r = [r.beta(1,1,n)*0.5 for n in n_loci]

    #one to use in the alternate approach (creates a lookup-array of the recombination values, essentially)
alt_r = samps = [[r.binomial(1,p,1000) for p in curr_r[i]] for i in range(len(curr_r))]


#vector of num of new gametes
#N = [10,100,1000]
N = [10,100,1000,10000]


def current_r_approach():
        return([r.binomial(1,curr_r[i]) for individ in range(n)])

def alternate_r_approach():
        return([r.choice(locus_array, size = n, replace = True) for locus_array in alt_r[i]])

def compare(func_1, func2):
    cur = {}

    alt = {}
    
    for n in N:
        cur[n] = []
        alt[n] = []
        print(n)
        for i in range(len(curr_gen_r)):
            print('\t%i'%i)
            #per suggestions at https://stackoverflow.com/questions/8220801/how-to-use-timeit-module, doing it the following way:
            REPEATS = 10
            NUMBER = 100
            cur[n].append(timeit.Timer(func_1).timeit(number=NUMBER)/NUMBER)
            #curr[n].append(min(timeit.Timer(current_approach).repeat(repeat=REPEATS, number=NUMBER))/NUMBER)
            alt[n].append(timeit.Timer(func_2).timeit(number=NUMBER)/NUMBER)
            #alt[n].append(min(timeit.Timer(alternate_approach).repeat(repeat=REPEATS, number=NUMBER))/NUMBER)

    return(cur, alt)



cur, alt = compare(current_r_approach, alternate_r_approach)



fig = plt.figure()
fig.suptitle('Time to create N new recombinant gametes of l loci, compared for\ncurrent (red) and proposed alternate (blue) recombination approaches')
for i,n in enumerate(N):
    ax = fig.add_subplot(2,2,i+1)
    ax.set_title('number of new recombinant gametes = %i' % n)
    ax.set_xlabel('number of recombining loci')
    ax.set_ylabel('time (sec)')
    ax.set_xscale("log", nonposx='clip')
    ax.plot(n_loci, cur[n], color = 'red')
    ax.plot(n_loci, alt[n], color = 'blue')

plt.show()







################################
# INDIVIDUALS QUERY RASTER CELLS
################################


#NOTE: NEVERMIND. No need to test this. Using %timeit on the command line the new approach is 
#clearly much faster: 1.33 ms per loop, vs 12.8 ms per loop, and this with only 500 individs and 3 scapes


def cur_query_approach():
    return([ind.query_habitat(land) for ind in self.individs.values()])

def alt_query_approach():
    return({i:[land.scapes[sc].raster[int(ind.y), int(ind.x)] for sc in land.scapes.keys()] for i,ind in pop.individs.items()})







