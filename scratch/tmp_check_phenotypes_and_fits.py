import numpy as np
import geonomics as gnx


def recalc_z(spp, trt, ind, monogen=False, burn=True):
    if burn:
        if not monogen:
            re_z = 0.5 + np.sum(np.mean(ind.g[trt.loci, :], axis=1) * trt.alpha)
        else:
            re_z = np.mean(ind.g[trt.loci, :], axis=1)[0]
    else:
        g = spp._get_genotypes(individs=[ind.idx])[ind.idx]
        if not monogen:
            re_z = 0.5 + np.sum(np.mean(g[trt.loci, :], axis=1) * trt.alpha)
        else:
            re_z = np.mean(g[trt.loci, :], axis=1)[0]
    return re_z

def check_z(spp, monogen=False, burn=True):
    for ind in spp.values():
        for trt_num, trt in spp.gen_arch.traits.items():
            z = ind.z[trt_num]
            re_z = recalc_z(spp, trt, ind, monogen=monogen, burn=burn)
            print(z, re_z)
            assert z == re_z
    return

def check_fit(spp, monogen=False, burn=True):
    for ind in spp.values():
        fit = ind.fit
        recalc_fit = 1
        for trt_num, trt in spp.gen_arch.traits.items():
            re_z = recalc_z(spp, trt, ind, monogen=monogen, burn=burn)
            e = ind.e[trt.lyr_num]
            recalc_trt_fit = 1 - (trt.phi * np.abs(e-re_z)**trt.gamma)
            recalc_fit *= recalc_trt_fit
        if len(spp.gen_arch.delet_loci) > 0:
            del_fit = gnx.ops.selection._calc_fitness_deleterious_mutations(spp,
                                                                  individs=[ind.idx])[0]
            recalc_fit *= del_fit
        print(fit, recalc_fit)
        assert fit == recalc_fit
    return

