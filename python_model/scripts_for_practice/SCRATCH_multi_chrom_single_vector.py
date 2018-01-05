
def get_chrom_breakpoints(l_c):
    breakpoints = np.array([0]+list(np.cumsum(sorted(l_c))[:-1]))
    assert np.alltrue(diff(bps) == np.array(sorted(l_c))), 'The breakpoints assigned will not produce chromosomes of the correct length'
    return(breakpoints)





def create_recomb_array(params):

    #get L (num of loci) and l_c (if provided; num of loci per chromsome) from params dict
    L = params['L']
    if ('l_c' in params.keys() and params['l_c'] <> None and len(params['l_c'] > 1)):
        l_c = params['l_c']
    
        #if l_c provided, check chrom lenghts sum to total number of loci
        assert sum(l_c) == L, 'The chromosome lengths provided do not sum to the number of loci provided.'

    #if params['recomb_array'] (i.e a linkage map) manually provided (will break if not a list, tuple, or np.array), 
    #then set that as the recomb_array, and check that len(recomb_array) == L
    if params['linkage_map'] <> None:
        recomb_array = np.array(params['recomb_array'])
        len(recomb_array) == L, "Length of recomb_array provided not equal to params['L']."

        return(recomb_array)

    
    #otherwise, create recomb array
    else:
        #start as an array of NaNs
        recomb_array = np.array([np.nan]*L)

        #set recomb values, using arguments provided in the params dict:

        #if a custom recomb_fn is provided, grab it
        if ('recomb_fn' in params['custom_fns'] and params['custom_fn']['recomb_fn'] <> None):
            recomb_fn = params['custom_fns']
            assert callable(recomb_fn), "The recomb_fn provided in params['custom_fns'] appear not to be defined properly as a callable function."
            #then call the draw_recomb_rate() function for each locus, using custom recomb_fn
            recomb_array = draw_recomb_rate(recomb_array, params, recomb_fn = recomb_fn)

        else:
            recomb_array = draw_recomb_rate(recomb_array, params) 


        #set recomb rate at the appropriate chrom breakpoints to 0.5
        bps = get_chrom_breakpoints(l_c)
        recomb_array[bps] = 0.5

        return(recomb_array)



