import numpy as np


def test_spp_get_genotypes(spp):
    if spp.gen_arch.use_tskit:
        res = [np.all(spp[idx].g== spp._get_genotypes(
            individs = [idx],
            biallelic=True)[0][spp.gen_arch.nonneut_loci,:]) for idx in [*spp]]
    else:
        res = [np.all(spp[idx].g == spp._get_genotypes(
            individs = [idx], biallelic=True)[0]) for idx in [*spp]]
    succeed = np.all(res)
    if succeed:
        return succeed
    else:
        return dict(zip([*spp], res))

#    idxs = [*spp]
#    gts = spp._get_genotypes()
#    both = dict(zip(idxs, gts))
#    checks = {}
#    for k, gt in both.items():
#        if not use_tskit:
#        #if not spp.gen_arch.use_tskit:
#            check = np.all(both[k] == spp[k].g)
#        else:
#            check = np.all(both[k][spp.gen_arch.nonneut_loci, :] == spp[k].g)
#        checks[k] = check
#    all_good = False not in [check == True for check in checks.values()]
#    if all_good:
#        return all_good
#    else:
#        print('UH OH')
#        return checks


def test_spp_get_phenotypes(spp):
    """
    Only works for 1 trait right now
    """
    if spp.gen_arch.use_tskit:
        man_z = [(0.5 + np.sum(np.mean(spp[idx].g, axis=1)
                         * spp.gen_arch.traits[0].alpha)) for idx in [*spp]]
    else:
        man_z = [(0.5 + np.sum(np.mean(spp[idx].g[spp.gen_arch.traits[0].loci, :], axis=1)
                         * spp.gen_arch.traits[0].alpha)) for idx in [*spp]]
    aut_z = [spp._get_z(individs=[idx]#, trait_num=0
                           )[0] for idx in [*spp]]
    succeed = np.allclose(man_z, aut_z)
    if succeed:
        return succeed
    else:
        return dict(zip([*spp], res))


