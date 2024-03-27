import geonomics as gnx
import tskit
import os

mod0 = gnx.make_model('GNX_tiny_model_params.py')
mod1 = gnx.make_model('GNX_tiny_model_params.py')
mods = [mod0, mod1]
tcs = []
for i, mod in enumerate(mods):
    reduced_n = [3, 5][i]
    mod.walk(10000, 'burn')
    mod.walk(50, 'main')
    spp = mod.comm[0]
    spp._reduce_to_size_n(reduced_n)
    print(spp)
    print('-'*80)
    tc = spp._tc
    tcs.append(tc)
    with open(f"mod{i}.svg", 'w') as f:
        f.write(tc.tree_sequence().draw_svg())
# add mod1 into mod0
tcs[0].union(tcs[1],
             # NOTE: no nodes in common in this instance,
             #       so just provide list of tskit.NULL
             #       equal in length to nodes in mod1's nodes btale
             node_mapping=[tskit.NULL]*tcs[1].nodes.num_rows,
             add_populations=True,
             check_shared_equality=True,
             record_provenance=True,
            )
with open('mod0_and_mod1.svg', 'w') as f:
    f.write(tcs[0].tree_sequence().draw_svg())
os.system('firefox mod0.svg mod1.svg mod0_and_mod1.svg')



