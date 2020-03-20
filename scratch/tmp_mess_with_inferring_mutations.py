import msprime
import tskit
import geonomics as gnx
import numpy as np

# make a simple gnx model
dir = '/home/drew/Desktop/stuff/berk/research/projects/sim/geonomics'
file = '/GNX_default_model_params.py'
mod = gnx.make_model(dir + file)
mod.walk(1000, 'burn')

# reduce the population to just 3
mod.comm[0]._reduce(10)

# randomly choose 3 loci and make them fixed for either 0 or 1
#fixed_loci = np.random.choice(range(10), 3, replace=False)
#print(fixed_loci)
#for ind in mod.comm[0].values():
#    ind.g[fixed_loci, :] = np.array([0,0,1,1,0,0]).reshape((3,2))
# NOTE: setting all individuals' genotypes to all 0s
for ind in mod.comm[0].values():
    ind.g = np.zeros((mod.comm[0].gen_arch.L, mod.comm[0].gen_arch.x))
new_p = mod.comm[0].gen_arch.p[:]
#new_p[fixed_loci] = [0, 1, 0]
mod.comm[0].gen_arch.p = new_p

# grab all the genotypes
 # (axes: 0 = individuals, 1 = loci, 2 = homologues)
genotypes = np.stack([ind.g for ind in mod.comm[0].values()])
#segregating_sites = np.where(genotypes.sum(axis=2).sum(
#                                    axis=0)/(2 * len(mod.comm[0])) % 1 != 0)[0]

# simulate a coalescent ancestry with number of samples equal to our species
ts = msprime.simulate(len(mod.comm[0]) * 2, Ne=1000,
                      length=mod.comm[0].gen_arch.L)

# mutate the simulated ancestry, with 3/gen mut rate and inf. sites
#ts = msprime.mutate(ts, rate=3, 
#                    model=msprime.InfiniteSites(msprime.NUCLEOTIDES),
#                    random_seed=2)

# then grab the tree
# NOTE: the TreeSequence object only has one, because no recombination was used
# in the sim
tree = ts.first()

# get the full TableCollection
tables = ts.dump_tables()

# set the sequence length
tables.sequence_length = mod.comm[0].gen_arch.L

# clear the mutations table
tables.mutations.clear()

# generate sufficient mutations, at only current nodes, to produce the starting
# 1-allele frequencies parameterized for a geonomics species
def make_mutations(spp, tables):
    print(tables.sites)
    print('ARE SITES ALL ANCESTRAL-STATE 0?')
    # get all the sample nodes, because those are where we want to plunk down
    # mutations
    #mutable_nodes = [n for n, flag in enumerate(
    #                                        tables.nodes.flags) if flag == 1]
    # get the starting frequencies for each site
    start_freqs = spp.gen_arch.p

    # get a set of all homologues, as tuples of (individual id, homologue idx)
    homologues = [*zip(np.repeat([*spp], 2),
                         [*range(spp.gen_arch.x)] * len(spp))]

    # make mutations for each site
    for site, freq in enumerate(start_freqs):
        # generate the number of mutations for this locus
        n_mutations = int(round(2 * len(spp) * freq, 0))
        print('n haploids')
        print(2*len(spp))
        print('n mutations')
        print(n_mutations)
        print('\n\n')
        # make sure we don't mutate either all or none, unless called for
        if n_mutations == len(spp) * 2 and freq < 1:
            n_mutations -= 1
        if n_mutations == 0 and freq > 0:
            n_mutations = 1
        #randomly choose and mutate n_mutations homologues from the population 
        curr_site_homologues = homologues[:]
        np.random.shuffle(curr_site_homologues)
        print('all homolgues')
        print(curr_site_homologues)
        homologues_to_mutate = curr_site_homologues[:n_mutations]
        print('homolgues to mutate')
        print(homologues_to_mutate)
        print('\n\n-------------')
        for ind, homol in homologues_to_mutate:
            spp[ind].g[site, homol] = 1
            node_id = spp[ind]._nodes_tab_ids[homol]
            tables.mutations.add_row(site, node=node_id, parent=-1,
                                     derived_state='1')
    return


site_ct = 0
# loop over all sites, so that each site's row goes into the sites table in
# order (such that 1.) there's no need to track how gnx sites
# map onto sites-table row ids, and 2.) there will be no need
# to deduplicate sites later on)
for site in range(mod.comm[0].gen_arch.L):
    #determine whether this is a neutral or non-neutral site
    if site in mod.comm[0].gen_arch.neut_loci:
        metadata='n'.encode('ascii')
    elif locus in mod.comm[0].gen_arch.nonneut_loci:
        metadata='t'.encode('ascii')
    # return a parsimonious set of state transitions, given the
    # observed genotypes and alleles for the samples in this tree
    #ancestral_state, mutations = tree.map_mutations(
    #            genotypes[:, site, :].flatten(), alleles = ('0', '1'))
    # add the variant's site to the sites table
    #tables.sites.add_row(position=site, ancestral_state=ancestral_state,
    tables.sites.add_row(position=site, ancestral_state='0',
                         metadata=metadata)
    # get the current n rows in the muts table
    parent_offset = len(tables.mutations)
    print('\n\nvar #', site_ct, 'parent offset', parent_offset)
    # add each mutation to the table
    #for i, mutation in enumerate(mutations):
        # grab the mutation's parent node
        #parent = mutation.parent
        # determine if this site has been subject to repeat mutations,
        # and if so then add the parent offset to the parent value, to
        # get the correct muts table id (i.e. row number) to be used
        # as the parent value
        #if parent != tskit.NULL:
        #    parent += parent_offset
        #print('\tsite #', site_ct, 'mut #', i+1)
        # add this mutation to the mutations table
        #tables.mutations.add_row(site, node=mutation.node,
        #                         parent=parent,
        #                         derived_state=mutation.derived_state)
                                 # NOTE: adding no metadata, to indicate
                                 # that these are fake mutations mapped to
                                 # the coalescent tree simulated by msprime
    print('..............................')
    site_ct +=1



# grab the nodes flags, which are 1 for current nodes, 0 for past nodes,
# into two separate objects
current_nodes = [*np.where(tables.nodes.flags == 1)[0]]
# reverse, so that I can pop nodes off the 'front'
current_nodes = current_nodes[::-1]
past_nodes = [*np.where(tables.nodes.flags != 1)[0]]
# create an empty list, to fill up the individual ids for each node
nodes_tab_individual_col = np.int32(np.ones(len(tables.nodes.flags))*-1)

# NOTE: there are no requirements or restrictions for the individuals and
# nodes tables (e.g. order, etc.), and thus the tables' order is not affected
# by the TableCollection.simplify algorithm.
# So, I could add either current or past individuals and nodes first;
# choosing to add past first, so that all 'real' individuals from the start
# of the geonomics simulation forward will wind up having their individuals and
# nodes rows in a single block at the tables' bottoms

# add an individual to the individuals table for
# each coalescent-simulated node before the current time
# NOTE: adding no metadata, and no location, to indicate that this is
# a 'fake' individual, invented just to match up to the nodes
# simulated for the starting population
for node in past_nodes:
    ind_id = tables.individuals.add_row(flags=0)
    # store its individual id in the nodes table's individuals column
    nodes_tab_individual_col[node] = ind_id

# create and add to the individuals table a new row for each real individual
for ind in mod.comm[0].values():
    # get the 'location' column info, which will include the x and y
    # positions of an individual, as well as the individual's
    # phenotypes and fitness, if traits are being used
    loc = [ind.x, ind.y]
    if mod.comm[0].gen_arch.traits is not None:
        loc = loc + ind.z + [ind.fit]
    # add a new row to the individuals table, setting the location
    # column's value to loc
    # NOTE: using the metadata column to store to the gnx
    # individual idx, for later matching to update
    # Individual._individuals_tab_id after tskit's simplify
    # algorithm filters individuals
    ind_id = tables.individuals.add_row(flags=1, location=loc,
        metadata=ind.idx.to_bytes(length=4, byteorder='little'))
    ind._individuals_tab_id = ind_id

    # assign the individual 2 randomly chosen nodes from the current time step,
    # and associate the individual's _individuals_tab_id with those 2 nodes
    # in some data structure to collect this
    ind_node_ids = [current_nodes.pop() for _ in range(mod.comm[0].gen_arch.x)]
    ind._set_nodes_tab_ids(*ind_node_ids)
    # add this individual's ind_id to the nodes_tab_individual_col list,
    # once for each node
    nodes_tab_individual_col[ind_node_ids] = ind_id
# make sure that all nodes were assigned to individuals
assert np.all(nodes_tab_individual_col >= 0), 'Some nodes not given individs'

# use that node-individual_ids data structure to reset the individuals column
# in the nodels table
nodes_cols = tables.nodes.asdict()
nodes_cols['individual'][:] = nodes_tab_individual_col
nodes_cols['time'] += 1
tables.nodes.set_columns(**nodes_cols)

# add mutations
make_mutations(mod.comm[0], tables)


# display the results
print(tables)
print('\n\n\n')
print(tables.tree_sequence().first().draw(format='unicode'))

mod.comm[0]._tc = tables


print('\nSTARTING FREQS:')
print(np.hstack([ind.g for ind in mod.comm[0].values()]).sum(axis=1)/(2*len(mod.comm[0])))
