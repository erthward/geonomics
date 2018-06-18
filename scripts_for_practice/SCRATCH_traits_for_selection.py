#scratch script used for trying to develop trait-level selection for polygenic traits

#additional packages
import operator

#get number traits and number 
n_trait = params['n_traits']
#NOTE: FIX THIS!!!!! I created a logical loop here, because I set the number of selective loci based on g.non_neutral, but then later act as though g.non_neutral has not yet been set, and set the values using the resulting s values...
n_selective = sum([ sum([int(x) for x in chrom]) for chrom in pop.genomic_arch.non_neutral.values()])



p_trait = 1/float(n_trait) 
p_trait = p_trait + r.normal(0, 0.35, n_trait)*p_trait - (p_trait**2)
diff = 1 - sum(p_trait)
p_trait = diff/n_trait + p_trait
trait_bins = list(r.multinomial(n_selective, p_trait ,size = 1)[0])
zeros = trait_bins.count(0) 
for zero in range(zeros):
    trait_bins[trait_bins.index(0)] = 1
    trait_bins[trait_bins.index(max(trait_bins))] = max(trait_bins) - 1

assert 0 not in trait_bins


#perhaps THEN SELECT TRAIT-LEVEL s VALUES FROM DIST (along lines of Thurman and Barret 2016 paper), THEN DISSOCIATE INTO COMPOSITE LOCUS-SPECIFIC s VALUES WITHIN EACH TRAIT, THEN RANDOMLY SCATTER THESE ACROSS AN OTHERWISE [0]*L array of s VALUES??

trait_bins = dict(zip(range(len(trait_bins)), trait_bins))
s = r.beta(0.5, 2, n_trait)
s_trait = dict(zip(range(n_trait), s))
s_locus = {}

for trait in s_trait.keys():
    s = list(r.dirichlet(0.5, trait_bins[trait])[0])
    chroms = []
    loci = []
    for s_val in s:
        success = False
        while not success:
            chrom = np.randint(0, params['n'])
            locus = np.randint(0, g.l_c[chrom])
            if g.non_neutral[chrom][locus] == False:
                success = True
                chroms.append(chrom)
                loci.append(locus)
                g.non_neutral[chrom][locus] = True
                g.s[chrom][locus] = s_val
                g.trait[chrom][locus] = trait
            else:
                pass


    s_locus[trait] = dict(zip(zip(chroms, loci), s_vals))




trait_order = [i[0] for i in sorted(s_trait.items(), key = operator.itemgetter(1), reverse = True)]






