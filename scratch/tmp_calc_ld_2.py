import numpy as np
import matplotlib.pyplot as plt

def tmp_calc_ld(speciome, plot = False):

    #TODO: 02-01-20: 
        # - debug what I wrote below to keep only seg sites
        # - aggain change scape size in the sweep params file to 20,20
        # - try this out and see if it fixes my issues
        # - if so, tweak linkage plot so that it labels axes by pos numbers for
        #     seg sites

    #TODO: I should also include (either as an alternative within this fn,
    #or as separate fn) the option to calculate D'

    #TODO: I keep getting warnings like the following, which could just be 
    #due to divison of small floating-point numbers, but I should figure out 
    #exactly what's going on and be sure everything checks out. WARNING:
    # stats.py:117: RuntimeWarning: invalid value encountered in double_scalars

    #speciome = _get_speciome(spp)
    n = np.shape(speciome)[0] #num individs
    x = np.shape(speciome)[2] #ploidy
    N = n*x
    #L = spp.gen_arch.L
    L = speciome.shape[1]
    assert L == np.shape(speciome)[1], ("The length of the 1st dimension "
                            "of speciome doesn't equal spp.genomic_arch.L")

    #get rid of fixed sites, which create problems for linkage calculation
    segregating = np.sum(speciome, axis = (1,2)) / (
                                    speciome.shape[1] * speciome.shape[2])
    segregating = np.int8([n != 0 and n != 1 for n in segregating])
    keep_inds = [i for i in range(L) if segregating[i] == 1]

    speciome = speciome[segregating, :, :]

    seg_L = len(segregating)

    r2_mat = np.zeros([seg_L]*2) * np.nan # set NaN as the "no data" value

    for i in range(seg_L):
        for j in range(i+1, seg_L):
            #calculate freq of allele 1 at locus i
            f1_i = np.sum(speciome[:,i,:], axis = None)/(N)
            #calculate freq of allele 1 at locus j
            f1_j = np.sum(speciome[:,j,:], axis = None)/(N)
            #calculate freq of chroms with 1_1 haplotype at loci i and j
            f11_ij = float(np.sum(speciome[:,[i,j],:].sum(axis = 1) ==2,
                                                        axis = None))/(N)
            D_1_1 = f11_ij - (f1_i * f1_j)
            r2 = (D_1_1**2)/(f1_i*(1-f1_i)*f1_j*(1-f1_j))
            if not 0 <= r2 <= 1:
                print("\t r2 was", r2)
                print("\t f1_i was", f1_i)
                print("\t f1_j was", f1_j)
                print("\t f11_ij was", f11_ij)
                print("\t D11 was ", D_1_1)
            else:
                print("r2 was", r2)
                print("f1_i was", f1_i)
                print("f1_j was", f1_j)
                print("f11_ij was", f11_ij)
                print("D11 was ", D_1_1)

            r2_mat[i,j] = r2
            r2_mat[j,i] = r2

    return(r2_mat)



def plot_linkage_vs_dist(r2_mat, recomb_rates, L=1001, sweep_loc=500):
    # TODO: Do I need to account for possibility of even number of recombination
    # events between loci when calculating distance? In other words, should
    # I convert the distance to a true cM distance, or just leave as the sum
    # of the interlocus recombination rates?
    r2s = []
    dists = []
    for i in range(L-1):
        for j in range(i, L-1):
            r2 = r2_mat[i, j]
            if not np.isnan(r2):
                dist = np.sum(recomb_rates[i+1:j+1])
                r2s.append(r2)
                dists.append(dist)
    fig = plt.figure()
    plt.plot(dists, r2s, '.r')
    plt.xlabel("recombination distance (sum of interlocus recomb. rates)")
    plt.ylabel("linkage ($R^{2}$)")
    plt.show()
    return(r2s, dists)
