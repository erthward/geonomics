import geonomics as gnx

import os
import shutil
import pandas as pd
import statsmodels.api as sm
import vcf

#set the data directory, and delete it if it already exists (so that we don't
#create mutliple, conflicting sets of data files)
data_dir = './GEONOMICS_mod-cline_params'
if os.path.isdir(data_dir):
    shutil.rmtree(data_dir)


mod = gnx.make_model('./test/validation/cline/cline_params/py')
#landscape and community will not be randomized between iterations, so I can
#just extract the non-neutral loci now
nonneutral_loci = mod.comm[0].gen_arch.traits[0].loci
mod.run(verbose = True)


#for each iteration
its_dirs = os.listdir(data_dir)
for it_dir in its_dirs:
    #read in the data
    files_dir = os.path.join(data_dir, it_dir, 'spp-spp_0')
    files = os.listdir(files_dir)
    vcf_file = [f for f in files if os.path.splitext(f)[1] == '.vcf']
    assert len(vcf_file) == 1
    vcf_file = vcf_file[0]
    vcf_reader = vcf.Reader(open(os.path.join(files_dir, vcf_file), 'r'))
    csv_file = [f for f in files if os.path.splitext(f)[1] == '.csv']
    assert len(csv_file) == 1
    csv_file = csv_file[0]
    csv = pd.read_csv(os.path.join(files_dir, csv_file))
    #grab all individuals' genotypes into an n_individs x n_loci array
    genotypes = np.ones((len(csv), spp.gen_arch.L))*99
    try:
        for loc in range(genotypes.shape[1]):
            rec = next(vcf_reader)
            for n, ind in enumerate(csv['idx']):
                genotype = sum([int(base) for base in rec.genotype(
                    str(ind)).data.GT.split('|')])/2
                genotypes[n, loc] = genotype
    except StopIteration:
        pass
    assert 99 not in genotypes
    glms = {}
    #parse the csv's 'e' column
    csv['env'] = [float(row.split(',')[0].lstrip('[')) for row in csv.e]
    #add each column of the genotypes array as a locus column in the csv
    for loc in range(genotypes.shape[1]):
        csv['loc%i' % loc] = genotypes[:,loc]
    #run a GLM for each locus
        try:
            glm = sm.GLM(csv['loc%i' % loc], csv['env'],
                family = sm.families.Binomial())
            glm_results = glm.fit()
            glms[loc] = glm_results
        except Exception as e:
            print(e)
    
    #grab all pvalues
    pvals = {loc:glm.pvalues[0] for loc, glm glms.items()} 
    loc = []
    pval = []
    for l, p, in pvals.items():
        loc.append(l)
        #divide by number of loci, as a simple correction for mutliple-testing
        pval.append(p/genotypes.shape[1])
    res = pd.DataFrame.from_dict({'loc':loc, 'pval':pval})
    #sort ascending, so that most significant loci should appear at top of res
    res = res.sort_values(by = 'pval')

#ANALYSIS IDEAS:

#try to pick out the correct (i.e. non-neutral) loci as outliers in a GLM
spp = mod.comm[0]
for locus in range(spp.gen_arch.L):
    logit_model = sm.GLM(gen_data[locus], env_data,
        family=sm.families.Binomial())


#also test hypothesis that mean fitness incrases and plateaus over model time
