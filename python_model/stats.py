#!/usr/bin/python
#stats.py

'''Functions for calculating common population-genetic statistics, and for doing so for specified stats at specified time
intervals during model run according params-dict contents.'''



#TODO: 
    # Turn the LD code below into functions: r2 and D' matrix functions
    # Move the het and maf functions that are embedded in population.birthday right now into here instead
    # Crosswalk this with arguments in params_dict (for which stats to calculate, and at what timsteps or freq)




from __future__ import division
from scipy.stats.stats import pearsonr
import time

populome = np.array([pop.individs[i].genome.genome[0] for i in pop.individs.keys()])
n = shape(populome)[0] #num individs
x = shape(populome)[2] #ploidy
N = n*x
L = pop.genomic_arch.L

r2_mat = np.zeros([L]*2)-1
alt_r2_mat = np.zeros([L]*2)-1

#first approach
start_1 = time.time()
for i in range(L):
    for j in range(i+1, L):
        loc_i_f1 = sum(populome[:,i,:])/(N)  #calculates freq of allele 1 at locus i
        loc_j_f1 = sum(populome[:,j,:])/(N)  #calculates freq of allele 1 at locus j
        loc_ij_f11 = sum(populome[:,[i,j],:].sum(axis = 1) ==2)/N #calculates freq of chroms with 1_1 haplotype at loci i and j
        D_1_1 = loc_ij_f11 - (loc_i_f1 * loc_j_f1)
        r2 = D_1_1**2/(loc_i_f1*(1-loc_i_f1)*loc_j_f1*(1-loc_j_f1))
        r2_mat[i,j] = r2
end_1 = time.time()



        
#alternate for calculating
start_2 = time.time()
for i in range(L):
    for j in range(i+1, L):
        g_i = list(populome[:,i,0])+list(populome[:,i,1]) #for locus i, concatenate list of second chroms to end of list of first chroms
        g_j = list(populome[:,j,0])+list(populome[:,j,1]) #for locus j, do the same
        alt_r2 = pearsonr(g_i, g_j)[0]**2
        alt_r2_mat[i,j] = alt_r2
end_2 = time.time()



print('FIRST APPROACH: %0.4f s' % (end_1-start_1))
print('SECOND APPROACH: %0.4f s' % (end_2-start_2))
print('SPEED APPROACH 2 IS %0.4fx SPEED APPROACH 1' % ((end_1-start_1)/(end_2-start_2)))


assert round(r2,5) == round(alt_r2, 5), 'ERROR: Appears the two estimations of r^2 differ at less than 5 sig digits precision.'




#plot of LD matrix
plt.imshow(np.clip(r2_mat, a_min = 0, a_max = None), interpolation = 'nearest')



#plot of mean linkage values
r2_list = [r2_mat[0,1]]
for i in range(1,L-1):
    r2_list.append(mean([r2_mat[i-1,i], r2_mat[i,i+1]]))
r2_list.append(r2_mat[L-2,L-1])

plt.scatter(range(L), r2_list, c = 'red', marker = 'o', s=25)

