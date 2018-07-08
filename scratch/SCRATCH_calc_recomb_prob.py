#!/usr/bin/python

import numpy.random as r

def p_recomb(R,i,j):
    seg = R[i:j]
    prob = [1-r for r in seg] 
    prob = 1 - product(prob) 
    return prob


#R = [0.1]*10
R = list(r.uniform(0,0.1,100))
p_0_50 = p_recomb(R,0,50)
p_50_100 = p_recomb(R,50,100)
p_0_100 = p_recomb(R,0,100)
assert round(p_0_100,10) == round(p_0_50 + p_50_100 - (p_0_50 * p_50_100),10), '\n\n\tThe values differ by %0.14f  ' % ((p_0_50 + p_50_100 - (p_0_50 * p_50_100)) - p_0_100)

#It works! This will always produce a proper probability. Yet notice that map distances can sum to much more than one along the chromosome:

print(sum(R))


def calc_map_length(R):
    return(sum(R))

def calc_avg_crossovers_per_meiosis(R):
    return(sum(R)
