#
#def create_recomb_lookup_array_and_index(genomic_arch, L = 10000):
#    #create lookup array:
#
#    #the recombination rates to be modeled (i.e. evenly spaced values from 0 to 0.5 inclusive, spaced at 0.0001 precision)
#    la_vals = np.array(list(np.arange(0,0.5,1e-4)) + [0.5])
#    #and a dictionary with the recomb rates as keys and their indices as values
#    la_dict = dict(zip(la_vals, range(len(la_vals))))
#    
#    #the number of 0/1 values to be used to model the recombination rate of each row (i.e. determines the
#    #precision to which a row can model a probability)
#    #now create the lookup array, with each row containing L 0 and 1 values whose proportional sum approximates r,
#    #i.e. the recombination rate modeled by that row
#    la = np.array([[0] * int(round(L-L*p,0)) + [1] * int(round(L*p,0)) for p in la_vals])
#    
#    
#    #use recombination vector to create a new recomb vector containing the lookup array row indices of the
#    #nearest 'modeled' recombination rate 
#    r_index = [la_dict[p] for p in [min(la_vals, key = lambda x: np.abs(x - val)) for val in genomic_arch.r]]
#    
#    return(la, r_index)
#
#
#
#
#
#def recombine(la, r_index, n_recombinants):
#
#    #draw recombination 'paths' for all new recombinants
#    #new_recombs = np.array([la[new_r[l],r.randint(0, L, num_recombs)] for l in range(len(new_r))])
#    new_recombs = np.array([r.choice(la[l,:], size = n_recombinants, replace = True) for l in r_index])
#
#    return(new_recombs)
#



########
# TAKE 2
########



#create a lookup array, for raster recombination of larger numbers of loci on the fly
    #NOTE: size determines the minimum distance between probabilities (i.e. recombination rates)
    #that can be modeled this way
def create_recomb_lookup_array(genomic_arch, size = 10000):

    
    #TODO: May need to work on this section below, to provide ability to create a list of lookup arrays, in case where number of loci being
    #modeled exceeds possible memory allotment for a single np.ndarray
    #then would loop over them by using int(locus/L) to increment along the list


    #determine how many arrays would be needed on this system for the given amount of loci
    #max_loci_per_array = 0
    #max_val = 100000
    #for i in [100, 1000, 10000, 100000]:
    #    try: 
    #        x = np.zeros([i,10000], dtype = np.int8)
    #    except MemoryError:
    #        max_val = i
    #        break
    #num_arrays = int(np.ceil(L/float(max_val))) 
    #arrays = [max_val] * num_arrays
    #arrays[-1] = pop.genomic_arch.L - max_val*(num_arrays)-1

    num_arrays = 1
    
   

    las = []
    for n in range(num_arrays):
        r_slice = pop.genomic_arch.r[n*max_val:(n*max_val + max_val)]
        la = np.zeros((len(r_slice),size), dtype = np.int8)

        for i, rate in enumerate(r_slice):
            la[i,0:int(round(size*rate))] = 1
        las.append(la)

    if len(las) >1:
        return(las)
    else:
        return(las[0])


def recombine(la, n_recombinants):
    return(np.array([r.choice(la[i,], size = n_recombinants, replace = True) for i in range(len(la))]))


