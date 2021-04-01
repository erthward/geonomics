import geonomics as gnx
import numpy as np


# FLAG DETERMINING WHETHER TO TEST TRAIT MUTATION OR DELETERIOUS MUTATION
mutate_trait = False

mod = gnx.make_model('./GNX_default_model_params.py')
mod.walk(10000, 'burn')
mod.walk(1)

spp = mod.comm[0]
ga = spp.gen_arch
re = ga.recombinations
trt = ga.traits[0]

off = [i.idx for i in spp.values() if i.age == 0] 

ga.mutables = [*ga.neut_loci]
np.random.shuffle(ga.mutables)

if mutate_trait:
    #PRINT STUFF BEFOREHAND
    print('ga.nonneut_loci', ga.nonneut_loci)
    print('trait loci', ga.traits[0].loci)
    print('trait locus index', ga.traits[0].loc_idx)
    print('unmutated genome:\n', spp[off[0]].g)
    print('mutated genome:\n', spp[off[-1]].g)
    nonneut_loci_b4 = set([*ga.nonneut_loci])


    gnx.ops.mutation._do_nonneutral_mutation(spp, [off[-1]], trait_nums=[0])
    #PRINT STUFF AFTERWARD
    print('trait loci', ga.traits[0].loci)
    print('trait locus index', ga.traits[0].loc_idx)
    print('unmutated genome:\n', spp[off[0]].g)
    print('mutated genome:\n', spp[off[-1]].g)
    nonneut_loci_af = set([*ga.nonneut_loci])
    new_locus = [*nonneut_loci_af.difference(nonneut_loci_b4)][0]

    mut_loc_idx = np.where(trt.loc_idx == np.where(
                                                trt.loci == new_locus)[0][0])[0][0]
    mut_homol = np.where(spp[off[-1]].g[mut_loc_idx, :] == 1)[0][0]
    #mutated_homologue = np.where(spp[off[-1]].g[np.where(
    #                        trt.loc_idx == np.where(
    #                            trt.loci == new_locus)[0][0])[0][0],:] == 1)[0][0]

    # MAKE SURE THE INFO IN GEONOMICS' NATIVE DATA STRUCTURES MATCHES THAT IN THE
    # TSKIT STRUCTURES
    print(spp._tc.mutations, '\n')
    print('last row should read:\n%i\t%i\t%i\t1\t-1' % (
                                    spp._tc.mutations.num_rows - 1,
                                    new_locus,
                                    spp[off[-1]]._nodes_tab_ids[ mut_homol]))

else:
    #PRINT STUFF BEFOREHAND
    print('delet loci', ga.delet_loci)
    print('unmutated genome:\n', spp[off[0]].g)
    print('mutated genome:\n', spp[off[-1]].g)
    nonneut_loci_b4 = set([*ga.nonneut_loci])

    gnx.ops.mutation._do_nonneutral_mutation(spp, [off[-1]], delet_s=0.1)

    #PRINT STUFF AFTERWARD
    print('delet loci', ga.delet_loci)
    print('unmutated genome:\n', spp[off[0]].g)
    print('mutated genome:\n', spp[off[-1]].g)
    nonneut_loci_af = set([*ga.nonneut_loci])
    new_locus = [*nonneut_loci_af.difference(nonneut_loci_b4)][0]

    mut_loc_idx = ga.delet_loc_idx[np.where(
                                    ga.delet_loci == new_locus)[0][0]]
    mut_homol = np.where(spp[off[-1]].g[mut_loc_idx, :] == 1)[0][0]
    #mutated_homologue = np.where(spp[off[-1]].g[np.where(
    #                        trt.loc_idx == np.where(
    #                            trt.loci == new_locus)[0][0])[0][0],:] == 1)[0][0]

    # MAKE SURE THE INFO IN GEONOMICS' NATIVE DATA STRUCTURES MATCHES THAT IN THE
    # TSKIT STRUCTURES
    print(spp._tc.mutations, '\n')
    print('last row should read:\n%i\t%i\t%i\t1\t-1' % (
                                    spp._tc.mutations.num_rows - 1,
                                    new_locus,
                                    spp[off[-1]]._nodes_tab_ids[ mut_homol]))


