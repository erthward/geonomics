import geonomics as gnx

import numpy as np
from copy import deepcopy
import shutil
import os
import vcf
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import statsmodels.api as sm

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

def calc_Fst(f0, f1, het0, het1, est_Hs = False):
    if f0 == f1:
        Fst = 0
    else:
        Ht = 2 * ((f0 + f1)/2) * (1- (f0+f1)/2)
        if est_Hs:
            Hs = (f0 * (1-f0)) + (f1 * (1-f1))
        else:
            Hs = np.mean([het0, het1])
        Fst = (Ht - Hs)/Ht
    return Fst

def test_island_model():
    #define island number, width, and diagonal distance
    n = 6
    w = 6
    d = 6
    #define the data directory, and remove it if it already exists (so that we
    #don't create multiple output data files in the same
    #directory and break the assertions below)
    data_dir = './GEONOMICS_mod-stepping_stone_params/'
    if os.path.isdir(data_dir):
        shutil.rmtree(data_dir)
    #create a Model from the params file
    print('\n\nMaking model...\n\n')
    mod = gnx.make_model(('./geonomics/tests/validation/stepping_stone/'
                          'stepping_stone_params.py'))
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
    #data structure to store islands' population sizes each timestep
    pop_sizes = {0:[], 1:[], 2:[], 3:[], 4:[], 5:[]}
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
        assert 0 not in old_vals, ("It appears "
            "some individuals landed in the 'sea' during the "
            "previous timestep but didn't die!")
        #record pop sizes
        for k, v in island_vals.items():
            pop_sizes[v].append(sum([old_val == k for old_val in old_vals]))
        old_island_labels = dict(zip([*mod.comm[0]],
            [island_vals[i] for i in old_vals]))
        mod.walk(mode = 'main', T = 1, verbose = True)
        #record the number of individuals whose island numbers
        #changed (i.e. who migrated this timestep)
        new_vals = [i[1] for i in mod.comm[0]._get_e()]
        try:
            assert 0 not in new_vals, ("It appears "
                "some individuals landed in the 'sea' during "
                "this timestep but didn't die!")
        except Exception as e:
            for k,v in mod.comm[0].items():
                print(v.fit)
                if v.e[1] == 0:
                    print('------')
                    print(k)
                    print(v.idx)
                    print(v.e)
                    print(v.x, v.y)
                    print(dict(zip([*mod.comm[0]], mod.comm[0].cells))[v.idx])
                    print(dict(zip([*mod.comm[0]], mod.comm[0].coords))[v.idx])
                assert k == v.idx, ('BIG PROBLEM: Not all individuals\''
                    'index values equal their dict keys.')
            print("PASSED THAT TEST")
            inds_and_new_vals = {i.idx: i.e[1] for i in mod.comm[0].values()}
            mod.plot(0,0)
            print("MADE inds_and_new_vals")
            inds = [k for k,v in inds_and_new_vals.items() if v == 0]
            mod.plot(0,0, individs = inds)
            print(inds)
            for k in inds:
                print(mod.comm[k].e)
            print(e)
            raise ValueError
        new_island_labels = dict(zip([*mod.comm[0]],
            [island_vals[i] for i in new_vals]))
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
    for subdir in os.listdir(data_dir):
        #data structures to hold Fst values at all timesteps
        all_Fst_var = {}
        all_Fst_HsHt = {}
        #identify unique timesteps for which data was saved
        timesteps = unique([int(re.search('(?<=t\-)\d+(?=_spp-spp)',
            f).group()) for f in os.listdir(os.path.join(
            data_dir, 'it--1/spp-spp_0'))])
        for step in timesteps:
            files = os.listdir(os.path.join(data_dir, subdir, 'spp-spp_0'))
            #read in VCF
            vcf_files = [f for f in files if os.path.splitext(f)[1] == '.vcf']
            vcf_file = [f for f in vcf_files if re.search(
                '(?<=t\-)%i(?=_spp)' % step, f)]
            assert len(vcf_file) == 1, ('More than one VCF!')
            vcf_file = vcf_file[0]
            vcf_reader = vcf.Reader(open(os.path.join(
                data_dir, subdir, 'spp-spp_0', vcf_file), 'r'))
            #read in CSV
            csv_files = [f for f in files if os.path.splitext(f)[1] == '.csv']
            csv_file = [f for f in csv_files if re.search(
                '(?<=t\-)%i(?=_spp)' % step, f)]
            assert len(csv_file) == 1, ('More than one CSV!')
            csv_file = csv_file[0]
            csv= pd.read_csv(os.path.join(data_dir, subdir, 'spp-spp_0', csv_file))
            # use individuals' second environmental value to pin them to their
            # islands
            csv['island'] = [island_vals[float(r[1].e.split(',')[1].rstrip(
                ']'))] for r in csv.iterrows()]
            #get list of individuals on each island
            island_lists = {i:[*csv[csv['island'] == i]['idx']] for i in range(n)}
            #get 1-allele frequencies and heterozygosities for each allele,
            #for each island
            freq = {}
            het = {}
            for loc in range(100):
                try:
                    rec = next(vcf_reader)
                    loc_dict = {}
                    het_dict = {}
                    for isl_num, isl_list in island_lists.items():
                        het_cts = 0
                        cts_allele_1 = []
                        for ind in isl_list:
                            ct_allele_1 = sum([int(base) for base in rec.genotype(
                                str(ind)).data.GT.split('|')])
                            cts_allele_1.append(ct_allele_1)
                            if ct_allele_1 == 1:
                                het_cts += 1
                        #divide the count of the 1 allele in each island pop by
                        #2 times the island pop's size (i.e. 2N)
                        pop_size = len(csv[csv['island'] == isl_num])
                        freq_allele_1 = sum(cts_allele_1)/(2*pop_size)
                        #and divide the count of heterozygotes by the island pop's
                        #size (i.e. N)
                        het_freq = het_cts/pop_size
                        loc_dict[isl_num] = freq_allele_1
                        het_dict[isl_num] = het_freq
                    freq[rec.POS] = loc_dict
                    het[rec.POS]= het_dict
                except Exception as e:
                    pass
            #use that info to calculate Fst between all pairs of islands
            #for each allele
            #and use mean migrations rates for all pairs of islands to calculate
            #expected Fst, by: Fst = 1/(4*N_e*m + 1)
            pop_pairs = set([tuple(sorted(list(i))) for i in migs.keys()])
            Fst_var = {}
            Fst_HsHt = {}
            exp_Fst = {}
            tot_mig_rates = {}
            for pop_pair in pop_pairs:
                Fst_var_list = []
                Fst_HsHt_list = []
                for loc in freq.keys():
                    pop_freqs = freq[loc]
                    pop_hets = het[loc]
                    f0 = pop_freqs[pop_pair[0]]
                    f1 = pop_freqs[pop_pair[1]]
                    het0 = pop_hets[pop_pair[0]]
                    het1 = pop_hets[pop_pair[1]]
                    mean_f = np.mean([f0, f1])
                    var_f = np.var([f0,f1])
                    #calculate Fst using var and means of allele freqs
                    #(see Hartl & Clark 2007, pg. 291)
                    Fst_var_val = var_f/(mean_f * (1-mean_f))
                    #if both pops have a freq of 0, division by 0 creates nan,
                    #so replace with Fst = 0
                    if np.isnan(Fst_var_val):
                        Fst_var_val = 0
                    Fst_var_list.append(Fst_var_val)
                    Fst_HsHt_val = calc_Fst(f0, f1, het0, het1, est_Hs = False)
                    Fst_HsHt_list.append(Fst_HsHt_val)
                mean_Fst_var = np.mean(Fst_var_list)
                mean_Fst_HsHt = np.mean(Fst_HsHt_list)
                Fst_var[pop_pair] = mean_Fst_var
                Fst_HsHt[pop_pair] = mean_Fst_HsHt
                mig_rates = [migs[pop_pair], migs[pop_pair[::-1]]]
                tot_mig_rate = np.sum(mig_rates)
                tot_mig_rates.update({pop_pair: tot_mig_rate})
                exp_Fst[pop_pair] = 1/((4*tot_mig_rate) + 1)

            #add Fst values to the data structures that hold Fst values for all
            #timesteps
            all_Fst_var[step] = Fst_var
            all_Fst_HsHt[step] = Fst_HsHt

            #plot the results, if this is the final timestep
            if step == max(timesteps):
                fig = plt.figure()
                plt.suptitle(("Validations test #2: stepping-stone model"
                    "\n%i timesteps") % mod.params.model.T)
                ax = fig.add_subplot(121)
                mod.plot(0,0)
                ax.set_title('Islands and their subpopulations')
                plt.xlabel('Landscape x')
                plt.ylabel('Landscape y')
                plt.title
                ax = fig.add_subplot(122)
                ax.set_title('IBD between island subpopulations')
                ax.set_xlim(0,6)
                x = []
                Fst_HsHt_y = []
                Fst_var_y = []
                exp_Fst_y = []
                tot_mig_rate_y = []
                for pop_pair in pop_pairs:
                    x.append(np.abs(pop_pair[0] - pop_pair[1]))
                    Fst_HsHt_y.append(Fst_HsHt[pop_pair])
                    Fst_var_y.append(Fst_var[pop_pair])
                    exp_Fst_y.append(exp_Fst[pop_pair])
                    tot_mig_rate_y.append(tot_mig_rates[pop_pair])
                plt.plot(x, Fst_HsHt_y, 'or')
                plt.plot(x, Fst_var_y, 'og')
                plt.plot(x, exp_Fst_y, 'ob')
                ax.set_xlabel('Pairwise interisland distance (scaled)')
                ax.set_ylabel('Pairwise genetic distance ($F_{ST}$)')
                ax.legend(labels=['$F_{ST}=\\frac{H_{T} - H_{S}}{H_{T}}$',
                                         '$F_{ST}=\\frac{var(p)}{mean(p)}$',
                                         '$E[F_{ST}]=\\frac{1}{1+4Nm}$'],
                                 loc = 'best',
                                 fontsize = 'medium')
                x = np.array(x).reshape((len(x), 1))
                Fst_HsHt_y = np.array(Fst_HsHt_y).reshape((len(Fst_HsHt_y), 1))
                Fst_var_y = np.array(Fst_var_y).reshape((len(Fst_var_y), 1))
                tot_mig_rate_y = np.array(tot_mig_rate_y).reshape(
                    (len(tot_mig_rate_y), 1))
                exp_Fst_y = np.array(exp_Fst_y).reshape((len(exp_Fst_y), 1))
                line_preds_x = np.linspace(1,5,1000).reshape((1000,1))
                #run and plot linear regressions for each, with alpha = 0.01
                marker_codes = ['-r', '-g', '-b']
                colors = ['red', 'green', 'blue']
                names = ['Fst=Ht-Hs/Ht', 'Fst=var(p)/mean(p)', 'Fst=1/(4Nm + 1)']
                for n, data in enumerate([Fst_HsHt_y, Fst_var_y, exp_Fst_y]):
                    #est = sm.OLS(data, x).fit()
                    est = sm.OLS(data, np.hstack((x, x**2))).fit()
                    line_preds = est.predict(np.hstack((line_preds_x,
                        line_preds_x**2)))
                    plt.plot(line_preds_x, line_preds, marker_codes[n])
                    r2 = est.rsquared
                    p = est.pvalues[0]
                    plt.text(line_preds_x[100], line_preds[100]+0.07, 'r=%0.3f' % r2,
                        color = colors[n])
                    plt.text(line_preds_x[100], line_preds[100]+0.05, 'p=%0.3f' % p,
                        color = colors[n])
                    plt.ylim(0,1.05)
                    #validate the results on basis of linear regression coefficient and
                    #p-value
                    assert p < 0.01, (("Regression test for %s at 0.01 significance "
                        "level.") % names[n])
                    assert est.params[0] > 0, (("Coefficient for %s was not "
                        "positive.") % names[n])
                #plot the inter-island migration rates on the same x axis
                ax2 = ax.twinx()
                plt.plot(x, tot_mig_rate_y, 'o', color = 'black')
                ax2.set_ylabel(('mean inter-island migration rate (individuals per '
                    'generation)'))
                plt.gca().legend(labels = ['mean inter-island migration rate'],
                    loc = 7, fontsize = 'medium')
                est = sm.OLS(np.log(tot_mig_rate_y), np.log(x)).fit()
                line_preds = est.predict(np.log(line_preds_x))
                plt.plot(line_preds_x, np.e**line_preds, '-', color = 'black')
                r2 = est.rsquared
                p = est.pvalues[0]
                plt.text(line_preds_x[100], line_preds[100]+0.07, 'r=%0.3f' % r2,
                    color = 'black')
                plt.text(line_preds_x[100], line_preds[100]+0.05, 'p=%0.3f' % p,
                    color = 'black')
                ylim((0, 1.1*max(tot_mig_rate_y)))
                plt.show()

    # plot pop sizes over time
    fig2 = plt.figure()
    plt.suptitle('Population size and Fst, over model time')
    ax1 = fig2.add_subplot(121)
    ax2.set_title('Population size, colored by island')
    ax2.set_xlabel('t')
    ax2.set_ylabel('pop size (N)')
    colors = ['red', 'orange', 'yellow', 'green', 'blue', 'purple']
    for k, v in pop_sizes.items():
        ax2.plot(range(len(v)), v, '-', color = colors[k])
    # plot Fst values for all timesteps
    ax2 = fig2.add_subplot(122)
    ax2.set_title('Fst, colored by inter-island distance')
    ax2.set_xlabel('t')
    ax2.set_ylabel('Fst')
    colors = ['#8baee5', '#4580e0', '#0056e2']
    cmap = LinearSegmentedColormap.from_list('my_cmap', colors, N=50)
    colors = cmap(np.linspace(0, 1, 5))
    x = sorted(all_Fst_HsHt.keys())
    max_Fsts = []
    for pop_pair in pop_pairs:
        y = [all_Fst_HsHt[t][pop_pair] for t in x]
        max_Fsts.append(max(y))
        ax2.plot(x, y, color=colors[np.abs(pop_pair[0] - pop_pair[1])-1],
                linewidth=3)
    ax2.set_ylim((0, 1.1 * max(max_Fsts)))
    ax2.set_xlim((0, 1.1*max(x)))
    ax2.legend(labels=[str(i) for i in range(1, n)], loc='best',
              title='Inter-island distance', fontsize='medium')
    plt.show()

    print('\nValidation test successful.\n')

    # return mod, migs, freq, Fst_var, Fst_HsHt, exp_Fst, N_est,
    # tot_mig_rate_y, x
    return

test_island_model()
