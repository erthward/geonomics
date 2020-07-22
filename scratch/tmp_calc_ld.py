import numpy as np

# flake8: noqa

def calc_ld(spp, spome, plot = False):

    #TODO: I should also include (either as an alternative within this fn,
    #or as separate fn) the option to calculate D'

    #TODO: I keep getting warnings like the following, which could just be 
    #due to divison of small floating-point numbers, but I should figure out 
    #exactly what's going on and be sure everything checks out. WARNING:
    # stats.py:117: RuntimeWarning: invalid value encountered in double_scalars

    speciome = spome
    n = np.shape(speciome)[0] #num individs
    x = np.shape(speciome)[2] #ploidy
    N = n*x
    L = spp.gen_arch.L
    assert L == np.shape(speciome)[1], ("The length of the 1st dimension "
                            "of speciome doesn't equal spp.genomic_arch.L")

    r2_mat = np.zeros([L]*2)-1 # -1 serves as a 'no data' value

    for i in range(L):
        for j in range(i+1, L):
            #calculate freq of allele 1 at locus i
            f1_i = np.sum(speciome[:,i,:], axis = None)/(N)
            #calculate freq of allele 1 at locus j
            f1_j = np.sum(speciome[:,j,:], axis = None)/(N)
            #calculate freq of chroms with 1_1 haplotype at loci i and j
            f11_ij = float(np.sum(speciome[:,[i,j],:].sum(axis = 1) ==2,
                                                        axis = None))/(N)
            D_1_1 = f11_ij - (f1_i * f1_j)
            r2 = (D_1_1**2)/(f1_i*(1-f1_i)*f1_j*(1-f1_j))
            if np.isnan(r2):
                r2 = 0
            r2_mat[i,j] = r2

    return(r2_mat)





def new_calc_ld(spome):
    L = spome.shape[1]
    mat = np.ones((L, L)) * np.nan
    for i in range(L):
        for j in range(i+1, L):
            print('now i = %i, j = %i' % (i, j))
            p1 = np.mean(spome[:, i, :])
            print('calcd p1')
            q1 = np.mean(spome[:, j, :])
            print('calcd q1')
            pq_11 = np.mean(spome[:, (i, j), :].sum(axis=1) == 2)
            print('calcd pq_11')
            D = pq_11 - (p1 * q1)
            print('calcd D')
            r2 = (D**2)/(p1 * (1 - p1) * q1 * (1 - q1))
            print('calcd r2')
            if (np.isnan(r2) or np.isinf(r2)):
                r2 = 0
                print('\tcoerced to 0')
            print('calcdr2')
            mat[i, j] = r2
            print('saved it')
            print('p1 = %0.4f' % p1)
            print('q1 = %0.4f' % q1)
            print('pq_11 = %0.4f' % pq_11)
            print('D = %0.4f' % D)
            print('r2 = %0.4f' % r2)
            print('----------------')
            assert 0 <= r2 <= 1, "r2 not in 0 - 1 interval!"
    return mat
