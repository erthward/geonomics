#!/usr/bin/python


import numpy as np
import numpy.random as r



#How to reassign r values to maintain crossover probabilities if locus B is inserted between A and C with
#prior recombination rate r_AC = 0.2

#e.g.      |----------------------------:-------------------------|
#          A           r_AB             B          r_BC           C




#choose new r_BC value such that r_AC prior to intervening mutation and afterward are the same
def calc_r_BC(r_AC, r_AB):
    return((r_AB-r_AC)/((2*r_AB)-1))


#not for use in model, just to demonstrate that the function above, based on algebra I worked out, performs
#correctly (i.e. preserves the prior recombination probability between loci A and C after inserting B and
#adjusting the two associated r values
def sim_test(r_AC):
    r_AB = r.uniform(0,r_AC,1)
    r_BC = calc_r_BC(r_AC, r_AB)
    n = 10000
    single_locus_prob = mean([sum(r.binomial(1,r_AC,n))/float(n) for i in range(1000)])
    two_locus_prob = mean([shape(np.where(np.array([r.binomial(1,R,n) for R in [r_AB, r_BC]]).sum(axis = 0) == 1))[1]/float(n) for i in range(1000)])
    assert round(single_locus_prob,3) == round(two_locus_prob,3)
    return(single_locus_prob, two_locus_prob)



#then would need to do the following:

#increment the trait loci for all traits at loci > mutated locus
def increment_trait_loci():
    #find all trait loci > mutated locus

    #add 1 to them

    #insert either 1 or 0 (depending whether mutation is neutral or not) into pop.genomic_arch.non_neutral

    pass



def parameterize_non_neutral_mutation():
    #for non-neutral mutations, use genome module functions and pop.genomic_arch.traits for the appropriate
    #trait to assign the mutation an alpha value (and possibly a dominance factor, if using that in the future?)

    pass


def insert_mutation():
    #insert mutation at appropriate locus in the mutated individual

    #insert 0s in all other individs

    #NOTE: see mutation.py for this, basically already done
    
    pass


