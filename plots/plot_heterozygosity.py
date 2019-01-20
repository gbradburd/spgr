import pyslim, msprime
import numpy as np
import spatial_slim as sps
from msprime import BranchLengthStatCalculator as bs
import os, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

usage = """
Usage:
    {} (script name)
""".format(sys.argv[0])

if len(sys.argv) != 2:
    raise ValueError(usage)

script = sys.argv[1]

def compute_heterozygosity(ts, time, num_targets, mutation_rate=1e-8):
    # mts = sps.SpatialSlimTreeSequence(msprime.mutate(ts, mutation_rate), dim=2)

    alive = ts.individuals_alive(time)
    targets = np.random.choice(np.where(alive)[0], num_targets, replace=False)
    target_node_list = []
    new_nodes = []
    k = 0
    for ind in targets:
        target_node_list.extend(ts.individual(ind).nodes)
        new_nodes.append([k, k+1])
        k += 2

    rts = ts.recapitate(1e-9, Ne=1000)
    sts = rts.simplify(target_node_list)
    bsc = msprime.BranchLengthStatCalculator(sts)
    new_targets = [sts.node(u[0]).individual for u in new_nodes]
    het = np.array([bsc.divergence([list(sts.individual(u).nodes)], [0.0, sts.sequence_length])[0][0]
            for u in new_targets])

    locs = ts.individual_locations()
    return (het, targets)


def plot_heterozygosity(ts, het, targets):
    locs = ts.individual_locations()
    scaled_het = (het - np.mean(het)) / np.std(het)
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    fig = plt.figure(figsize=(6, 6 * ymax / xmax))
    colors = ['c' if h > 0 else 'm' for h in scaled_het]
    ax = fig.add_subplot(111)
    ax.scatter(locs[targets, 0], locs[targets, 1],
               s=400 * np.abs(scaled_het),
               alpha=0.75,
               c=colors)
    return fig


# for script in ("diff_valleys.slim", "valleys.slim", "flat_map.slim"):
# # script should be on the command line

num_gens = 10
treefile = sps.run_slim(script = script,
                        seed = 23, 
                        SIGMA = 1.0,
                        W = 50.0, 
                        K = 5.0,
                        NUMGENS = num_gens,
                        BURNIN=10000)
outbase = ".".join(treefile.split(".")[:-1])

ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

num_targets = 100
het_time = 0
hetfile = outbase + ".heterozygosity.txt"
if os.path.isfile(hetfile):
    print(hetfile, "already exists.")
    data = np.loadtxt(hetfile)
    het = data[:, 0]
    targets = np.array([np.int(u) for u in data[:, 1]])
    assert(max(np.abs(targets - data[:, 1])) == 0)
    if len(targets) != num_targets:
        print("Number of targets does not match saved file" + hetfile)
else:
    print(hetfile, "does not exist, computing.")
    het, targets = compute_heterozygosity(ts, het_time, num_targets)
    data = np.column_stack([het, targets])
    np.savetxt(hetfile, data)

fig = plot_heterozygosity(ts, het, targets)
fig.savefig(outbase + ".heterozygosity.pdf")
plt.close(fig)



