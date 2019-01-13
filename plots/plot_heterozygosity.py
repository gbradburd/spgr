import pyslim, msprime
import numpy as np
import spatial_slim as sps
from msprime import BranchLengthStatCalculator as bs

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_heterozygosity(ts, time, num_targets, mutation_rate=1e-8):
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
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    fig = plt.figure(figsize=(6, 6 * ymax / xmax))
    colors = ['b', 'c', 'm', 'y']
    ax = fig.add_subplot(111)
    ax.scatter(locs[targets, 0], locs[targets, 1],
               s=het/4,
               alpha=0.75,
               c='black')
    return fig


for script in ("valleys.slim", "flat_map.slim"):
    num_gens = 10
    treefile = sps.run_slim(script = script,
                            seed = 23, 
                            SIGMA = 1.0,
                            W = 50.0, 
                            K = 5.0,
                            NUMGENS = num_gens,
                            BURNIN=100)
    outbase = ".".join(treefile.split(".")[:-1])

    ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

    fig = plot_heterozygosity(ts, 0, 100)
    fig.savefig(outbase + ".heterozygosity.pdf")
    plt.close(fig)



