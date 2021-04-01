import geonomics as gnx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tskit

# get a model
fig = plt.figure()
mod = gnx.run_default_model()
fig.show()
mod.walk(2500)

# get the model's sorted TableCollection
spp = mod.comm[0]
tc = spp._tc
tc.sort()

# grab the edges and nodes as DataFrames
edgedf = pd.DataFrame(tc.edges.asdict())
nodedf = pd.DataFrame({k:v for k, v in tc.nodes.asdict(
            ).items() if k in ['time', 'individual']})

# and plot a histogram of the total branch lengths
# for all the trees in the TreeSequence
fig2 = plt.figure()
ax1 = fig2.add_subplot(311)
ts = tc.tree_sequence()
ts_it = ts.trees()
lens = []
for t in ts_it:
    lens.append(t.get_total_branch_length())
plt.hist(lens, bins=50)

# tables identical after simplification, if no individuals marked as samples
tc.simplify()
edgedf2 = pd.DataFrame(tc.edges.asdict())
nodedf2 = pd.DataFrame({k:v for k, v in tc.nodes.asdict(
            ).items() if k in ['time', 'individual']})
print('edges same?', np.all(edgedf == edgedf2))
print('nodes same?', np.all(nodedf == nodedf2))

# draw first tree
ts = tc.tree_sequence()
t = ts.first()
print(t.draw(format='unicode'))

# compare total branch lengths
ax2 = fig2.add_subplot(312)
ts = tc.tree_sequence()
ts_it = ts.trees()
lens = []
for t in ts_it:
    lens.append(t.get_total_branch_length())
plt.hist(lens, bins=50)


# get all current individuals and mark them as flag == 1 (i.e. set them as
# samples)
node_ids = np.hstack([[*ind._node_ids.values()] for ind in spp.values()])
#flags = tc.nodes.flags
#for s in samp:
#    flags[s] = tskit.NODE_IS_SAMPLE
#tc.nodes.set_columns(flags=flags, population=tc.nodes.population,
#                     time=tc.nodes.time)

# now simplify again
new_node_ids = tc.simplify(node_ids)

#tables are still equal to what they were
edgedf3 = pd.DataFrame(tc.edges.asdict())
nodedf3 = pd.DataFrame({k:v for k, v in tc.nodes.asdict(
            ).items() if k in ['time', 'individual']})
try:
    print('edges same?', np.all(edgedf == edgedf3))
    print('nodes same?', np.all(nodedf == nodedf3))
except ValueError as e:
    print('Data frames are different!\n\nHere are the original node_ids: %s\n\n And here are the new ones: %s\n\n' % (node_ids, new_node_ids))

# but again check the total branch lengths
ax3 = fig2.add_subplot(313)
ts = tc.tree_sequence()
ts_it = ts.trees()
lens = []
for t in ts_it:
    lens.append(t.get_total_branch_length())
plt.hist(lens, bins=50)

fig2.show()

# draw first tree again
ts = tc.tree_sequence()
t = ts.first()
print(t.draw(format='unicode'))

