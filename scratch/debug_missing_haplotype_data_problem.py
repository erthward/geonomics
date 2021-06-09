import pandas as pd

#hap dict
hap_dict = {n:hap for n,hap in enumerate([*ts.haplotypes()])}

# inds dict (by ordinal number and idx)
ind_dict = {n:ind for n, ind in enumerate([*spp])}

# get first tree is tree seq (i.e. tree for locus 0)
spp = mod.comm[0]
spp._sort_simplify_table_collection()
tc = spp._tc
ts = tc.tree_sequence()
tree = next(ts.trees())
# affirm that it's locus 0
assert tree.index == 0

# get samples with missing data (i.e. isolated samples) for this locus
missing = [u for u in tree.samples() if tree.is_isolated(u)]

# get the ids of the individs with missing data for this locus
individ_idxs_w_missing = [ind_dict[int(n/2)] for n in missing]
individs_w_missing = [spp[idx] for idx in individ_idxs_w_missing]

# get the first individual with missing data for this locus
ind = individs_w_missing[2]

# get the nodes table ids for this individ
nodes_ids = [*ind._nodes_tab_ids.values()]


def get_edges_df(tc):
    """
    function to get pandas df from first 4 cols of edges table
    """
    l = tc.edges.left
    r = tc.edges.right
    p = tc.edges.parent
    c = tc.edges.child
    assert len(l) == len(r) == len(p) == len(c)
    df = pd.DataFrame.from_dict({'left': l,
                                 'right': r,
                                 'parent': p,
                                 'child': c})
    return df


# get a pd df of the edges table's first 4 cols
edf = get_edges_df(tc)

# get the rows of the edges df that pertain to this individ's nodes
# and to this locus
edf_this_ind = edf.iloc[(np.int8(edf.child == nodes_ids[0]) +
                         np.int8(edf.child == nodes_ids[1])) == 1, :]
