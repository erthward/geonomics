import numpy as np
import time
import bisect

def m1(sh, event, locs, L):
    out = np.array([(sh + bisect.bisect_left(event,
                                             loc) % 2) % 2 for loc in locs])
    return out


# ~1.5x slower
def m2(sh, event, locs, L):
    #NOTE: doesn't yet work where locus and recomb breakpoint are identical,
    # but that doesn't matter because it's already 1.5x slower
    aug_event = [0] + [*event] + [L-1]
    out = np.hstack([np.repeat(i%2,
                               aug_event[i+1] - aug_event[i]) for i in range(
                                                    len(aug_event) - 1)])[locs]
    return out


# ~ 2x slower
def m3(sh, event, locs, L):
    aug_event = np.array([*event] + [L-1])
    out = np.array([(aug_event >= loc).argmax() % 2 for loc in locs])
    if sh:
        out = (out + 1) % 2
    return out




##########################
# I think I need a memoized version 

# here's a go...

# build the cache
def make_cache(L, locs):
    cache = dict(zip(range(L), [sum(np.array(locs) <= l) for l in range(L)]))
    return cache

# update the cache after a mutation
def update_cache(cache):
    pass

def m4(sh, event, locs, L, cache):
    if len(event) == 0:
        subsetter = [sh]*len(locs)
    else:
        pre_locs = np.array([0] + [cache[e] for e in event] + [cache[L-1]])
        diffs = pre_locs[1:] - pre_locs[:-1]
        subsetter = np.int8(np.hstack([[i % 2] * n for i, n in enumerate(
                                                                      diffs)]))
        if sh:
            subsetter = [(n+1) % 2 for n in subsetter]
    return(subsetter)




############################
# actually may need to just cache all recombination events beforehand,
# using bitarray subsetters like before

