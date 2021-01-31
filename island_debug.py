import vcf
import numpy as np
import pandas as pd


filename = "./GNX_mod-island_params/it--1/spp-spp_0/mod-island_params_it--1_t-770_spp-spp_0"
n_inds = None


def get_real_gt(idx):
    nextline = False
    gt_list = []
    with open('%s.fasta' % filename, 'r') as f:
        for line in f:
            if not nextline and line.startswith('>'+str(idx)):
                nextline = True
            elif nextline:
                gt_list.append(line)
                nextline=False
            else:
                pass
    gt_array = np.stack([[int(n) for n in hap if n != '\n'] for hap in gt_list]).T
    return gt_array


def get_rec_gt(idx):
    vr = vcf.Reader(open('%s.vcf' % filename, 'r'))
    rec_gt_arr = np.nan * np.ones((100,2))
    for loc in range(100):
        rec = next(vr)
        gt = [int(n) for n in rec.genotype(str(idx)).data.GT.split('|')]
        rec_gt_arr[loc, :] = np.int8(gt)
    return rec_gt_arr

df = pd.read_csv('%s.csv' % filename)
idxs = df.idx.values
if n_inds is None:
    n_inds = len(idxs)
run_idxs = np.random.choice(idxs, size=n_inds, replace=False)

res_dict = {}
print("COMPARING FOR %i INDIVIDUALS" % n_inds)
for idx in run_idxs:
    print('IDX: %i' % idx)
    real_gt = get_real_gt(idx)
    rec_gt = get_rec_gt(idx)
    res_dict[idx] = np.all(real_gt == rec_gt)
