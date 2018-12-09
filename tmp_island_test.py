import numpy as np
from copy import deepcopy
import geonomics as gnx
import os
import vcf
import pandas as pd

def make_island_landscape(n, w, d):
    """
    n = number of islands
    w = island width
    d = interisland diagonal distance (i.e. distance/sqrt(2))
    """
    #determine land dimension (i.e. nrows & ncols, hence n diagonal cells)
    dim = int(n * (w + d) + d)
    #create landscape
    land = np.zeros([dim]*2)
    island_ul_corners = [int(i) for i in np.linspace(d, dim-w-d, n)]
    for ulc in island_ul_corners:
        land[ulc:ulc+w, ulc:ulc+w] = 1
    #create a second raster, where each island's hab values are its island
    #number (to be used later to estimate migration rate, to then compare
    #results to expected results of stepping-stone model
    island_labels = deepcopy(land)
    for ulc in island_ul_corners:
        island_labels[ulc:ulc+w, ulc:ulc+w] = (ulc/(w+d)+1)/n

    return land, island_labels

def test_island_model():
    #define island number, width, and diagonal distance
    n = 6
    w = 6
    d = 6
    #create a Model from the params file
    print('\n\nMaking model...\n\n')
    mod = gnx.make_model('./test/validation/island/island_params.py')
    #replace Layer 0's raster with an island raster (6 islands, each 6
    #cells wide and 6 cells diagonally spaced), and Layer 1's raster
    #with a raster containing island labels
    islands, island_label_rast = make_island_landscape(n, w, d)
    mod.land[0].rast = islands
    mod.land[1].rast = island_label_rast
    #reset the Species' carrying capacity raster
    mod._set_K(0, mod.land)
    mod.comm[0]._do_movement(mod.land)
    #create a lookup-dict for the island numbers of each non-zero value in
    #island_labels
    island_vals =  dict(zip(sorted([*np.unique(island_label_rast)])[1:],
        range(n)))
    #run the model
    #burn the model in
    mod.walk(mode = 'burn', T = 1000000, verbose = True)
    #empty dict to hold proportions of each timestep's population
    #that underwent each of the possible inter-island pairwise
    #migration events (i.e. island 0 to 1; 1 to 0; 0 to 2; etc...)
    migs = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                migs.update({(i,j):[]})
    #walk Model timestep by timestep
    for t in range(mod.params.model.T):
        #create a dictionary to keep count of how many individuals make
        #each of the possible migration events
        mig_counts = {}
        for i in range(n):
            for j in range(n):
                if i != j:
                    mig_counts.update({(i,j):0})
        #record individuals' current island locations
        #(i.e. their second habitat values)
        old_vals = [i[1] for i in mod.comm[0]._get_e()]
        old_island_labels = dict(zip([*mod.comm[0]],
            [island_vals[i] for i in old_vals]))
        assert 0 not in old_vals, ("It appears "
            "some individuals landed in the 'sea' during the "
            "previous timestep and didn't die!")
        mod.walk(mode = 'main', T = 1, verbose = True)
        #record the number of individuals whose island numbers
        #changed (i.e. who migrated this timestep)
        new_vals = [i[1] for i in mod.comm[0]._get_e()]
        new_island_labels = dict(zip([*mod.comm[0]],
            [island_vals[i] for i in new_vals]))
        assert 0 not in new_vals, ("It appears "
            "some individuals landed in the 'sea' during "
            "this timestep and didn't die!")
        not_newborns = [ind for ind in [
            *old_island_labels] if ind in [*new_island_labels]]
        #get number of individuals who underwent each of the possible
        #migration events 
        for ind in not_newborns:
            if new_island_labels[ind] != old_island_labels[ind]:
                mig_event = (old_island_labels[ind], new_island_labels[ind])
                mig_counts[mig_event] += 1
        #append the proportion of individuals who underwent each
        #migration event to the mig_props dict
        [migs[mig_event].append(
            ct) for mig_event, ct in mig_counts.items()]
            #ct/len(not_newborns)) for mig_event, ct in mig_counts.items()]
    #get mean migration rates
    migs = {k:np.mean(v) for k,v in migs.items()}
    #divide mean population size by 6 to get estimated per-island pop sizes
    #NOTE: This needn't necessarily be close to N_e...
    N_est = np.mean(mod.comm[0].Nt)/n

    #for each subdirectory
    data_dir = './GEONOMICS_mod-island_params/'
    for subdir in os.listdir(data_dir):
        files = os.listdir(os.path.join(data_dir, subdir, 'spp-spp_0'))
        #read in VCF
        vcf_file = [f for f in files if os.path.splitext(f)[1] == '.vcf']
        print(vcf_file)
        assert len(vcf_file) == 1
        vcf_file = vcf_file[0]
        vcf_reader = vcf.Reader(open(os.path.join(
            data_dir, subdir, 'spp-spp_0', vcf_file), 'r'))
        #read in CSV
        csv_file = [f for f in files if os.path.splitext(f)[1] == '.csv']
        print(csv_file)
        assert len(csv_file) == 1
        csv_file = csv_file[0]
        csv= pd.read_csv(os.path.join(data_dir, subdir, 'spp-spp_0', csv_file))
        # use individuals' second environmental value to pin them to their
        # islands
        csv['island'] = [island_vals[float(r[1].e.split(',')[1].rstrip(
            ']'))] for r in csv.iterrows()]
        #get list of individuals on each island
        island_lists = {i:[*csv[csv['island'] == i]['idx']] for i in range(n)}
        #get 1-allele frequencies for each allele, for each island
        freq = {}
        for loc in range(100):
            try:
                rec = next(vcf_reader)
                loc_dict = {}
                for isl_num, isl_list in island_lists.items():
                    cts_allele_1 = []
                    for ind in isl_list:
                        ct_allele_1 = sum([int(i) for i in rec.genotype(
                            str(ind)).data.GT.split('|')])
                        cts_allele_1.append(ct_allele_1)
                    #divide the count of the 1 allele in each island pop by
                    #2 times the island pop's size (i.e. 2N)
                    freq_allele_1 = sum(cts_allele_1)/(2*len(
                        csv[csv['island'] == isl_num]))
                    loc_dict[isl_num] = freq_allele_1
                freq[rec.POS] = loc_dict
            except Exception as e:
                pass
        #use that info to calculate Fst between all pairs of islands
        #for each allele
        #and use mean migrations rates for all pairs of islands to calculate
        #expected Fst, by: Fst = 1/(4*N_e*m + 1)
        pop_pairs = set([tuple(sorted(list(i))) for i in migs.keys()])
        Fst = {}
        exp_Fst = {}
        for pop_pair in pop_pairs:
            Fst_list = []
            for loc in freq.keys():
                pop_freqs = freq[loc]
                f0 = pop_freqs[pop_pair[0]]
                f1 = pop_freqs[pop_pair[1]]
                mean_f = np.mean([f0, f1])
                var_f = np.var([f0,f1])
                Fst_val = var_f/mean_f
                #if both pops have a freq of 0, division by 0 creates nan,
                #so replace with Fst = 0
                if np.isnan(Fst_val):
                    Fst_val = 0
                Fst_list.append(Fst_val)
            mean_Fst = np.mean(Fst_list)
            Fst[pop_pair] = mean_Fst
            mig_rates = [migs[pop_pair], migs[pop_pair[::-1]]]
            tot_mig_rate = np.sum(mig_rates)
            exp_Fst[pop_pair] = 1/((4*tot_mig_rate) + 1)

    return mod, migs, freq, Fst, exp_Fst, N_est

