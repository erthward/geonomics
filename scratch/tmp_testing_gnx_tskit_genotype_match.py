import numpy as np
import tskit
import geonomics as gnx

def check_genotypes_match(spp):
    tc = spp._tc
    try:
        ts = tc.tree_sequence()
    except Exception:
        tc.sort()
        ts = tc.tree_sequence()
    samps = []
    for ind in spp.values():
        samps.append(ind._nodes_tab_ids[0])
        samps.append(ind._nodes_tab_ids[1])
    gm_genotypes = []
    gn_genotypes = []
    for site, var in enumerate(ts.variants(samples=samps)):
        gm_gt_to_append = var.genotypes
        gm_genotypes.append(gm_gt_to_append)
        gn_gt_to_append = np.hstack([ind.g[site, :] for ind in spp.values()])
        gn_genotypes.append(gn_gt_to_append)
        #print(gm_gt_to_append == gn_gt_to_append)
    gm = np.vstack(gm_genotypes)
    gn = np.vstack(gn_genotypes)
    #print('GN')
    #print(gn.shape)
    #print(gn)
    #print('GM')
    #print(gm.shape)
    #print(gm)
    return(np.all(gm == gn))


def check_if_sites_fixed(spp):
    speciome = np.hstack([ind.g for ind in spp.values()])
    freq = speciome.sum(axis=1)/(2*len(spp))
    print(np.hstack([ind.g for ind in spp.values()]).sum(axis=1)/(2*len(spp)))
    return(0 in freq)


def run_checks(mod, n, until_fixed=False):
    spp = mod.comm[0]
    match_res = []
    fix_res = []

    if not until_fixed:
        for _ in range(n):
            print('-----------------\n')
            mod.walk(1)
            match = check_genotypes_match(spp)
            fixed = check_if_sites_fixed(spp)
            print('MATCH:', match)
            print('FIXED:', fixed)
            match_res.append(match)
            fix_res.append(fixed)
    else:
        not_fixed = True
        while not_fixed:
            print('-----------------\n')
            mod.walk(1)
            match = check_genotypes_match(spp)
            fixed = check_if_sites_fixed(spp)
            print('MATCH:', match)
            print('FIXED:', fixed)
            match_res.append(match)
            fix_res.append(fixed)
            not_fixed = not fixed
    return(match_res, fix_res)
