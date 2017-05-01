#!/usr/bin/python

import numpy as np






#Get genotype
#NOTE: Incorporate h here, so that dominance, overdominance, codominance, etcetera can be modeled...






#Get the relative viabilities (proportional probabilities of survival) for a given environmental value and a
#list of genotypes at a certain locus for the individuals found there
def get_viability(e, g):
    #the absolute value of the difference between the env value and the genotype
    return(np.abs(e-np.array(g)))

#Get the vector of mortalities (probabilies of death) for a given density-dependent Pr(death) at a cell, the
#environmental value at that cell, a list of genotypes at a certain locus for the individuals found there, and the selection
#coefficient at that locus
    #This is calculated as: SIGMA_n->N,_ij d_ij*(1+diff),
    #where N is the number of individuals in cell i,j, d is the background density-dependent Pr(death) for
    #that cell, and diff is a correction factor, centered on a 0.5 difference between genotype and environ and
    #scaled so that the ratio of min to max possibile viabilities at this cell == 1-s, per common diploid
    #popgen theory

def get_prob_death(d,e,g,s):
    #print('e')
    #print(e)
    #print('d')
    #print(d)
    #print('g')
    #print(g)
    #print('s')
    #print(s)
    #c is the correction factor used to ensure that the resulting mortalities 
    c = -1*s/(2-s)
    #print('c')
    #print(c)
    #get the relative viabilities
    w = get_viability(e, g)
    #print('w')
    #print(w)

    #NOTE: IMPORTANT QUESTION: Should competition be completely relative (e.g. 4 maximally unfit individuals
    #in a cell would have equal chances of survival as four equally moderately fit individuals or 4 maximally
    #fit individuals; the way I have it coded now) or would those individuals perform worse even in the absence
    #of better-fit conspecifics in the same cell?? Not actually sure HOW important this is, because might only
    #be a matter of time before those individuals move into a neighboring cell, or have neighbors move in, and
    #thus the competitive environment changes...

    #get the scaled differential that will be used to calculate mortalities
    diff = ((0.5-w)/0.5)*c
    #print('diff')
    #print(diff)
    #calculate mortalities
    d_ind = d*(1+diff)
    #print('mort')
    #print(mort)
    
    #Check that the ratio of the min to max mortalities equals 1 - (ratio of max to min genotype values)*s
    assert np.allclose(d_ind.min()/d_ind.max(), 1-((max(g) - min(g))*s))

    #Now, tweak it if need be to ensure that it stays between 0 and 1
    #NOTE: NEED TO CONSIDER IF THERE IS A BETTER, LESS BIASING WAY OF DOING THIS
    d_ind = np.array([max(min(0.999999, val), 0.000001) for val in d_ind])
    print('mort again')
    print(d_ind)



    #Check that the mean of the mortalities equals the background density-dependent Pr(death) for the given cell
        #NOTE: at s = 0, the mean is precisely equal to the background d, but it diverges toward s = 1
        #I think this is a result of numerical computational error of margin, but need to double check (and if
        #it is, not yet sure how to rectify it)
    #NOTE: Actually, leaving this out for now, because of course the mean won't always equal d since the
    #vector w won't always be balanced. But for the moment operating on the assumption that cells where the
    #mean overshoots and cells where it undershoots will balance out on the average, and are a more realistic
    #reflection of fact that fitness is both competitive (i.e. relative, i.e. fitter than your neighbor) and
    #performative (i.e. absolute, i.e. fit enough to do well in your place, regardless of your neighbor)
    #assert np.allclose(np.mean(mort), d)

    
    
    #NOTE: NEED TO FIURE OUT BEST WAY TO KEEP MORTALITIES BETWEEN 0 AND 1
    #Check that 0<= mort <= 1
    assert np.alltrue(d_ind >= 0)
    assert np.alltrue(d_ind <= 1)

    #return
    return(d_ind)

    



