import pyslim, msprime
import numpy as np
import spatial_slim as sps
from msprime import BranchLengthStatCalculator as bs

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_isolation_by_distance(ts, targets, num_indivs):
    """
    Compute and plot mean sequence divergence of num_individuals to each of num_targets.
    """
    # a snapshot of the individual locations
    locs = ts.individual_locations()
    alive = ts.individuals_alive(0)
    indivs = list(np.random.choice(np.where(alive)[0], num_indivs, replace=False))
    all_indivs = targets + indivs
    num_targets = len(targets)

    # simplify so this goes faster
    individual_node_list = []
    simplified_node_list = []
    k = 0
    for ind in all_indivs:
        individual_node_list.extend(ts.individual(ind).nodes)
        simplified_node_list.append([k, k+1])
        k = k + 2

    sts = ts.simplify(individual_node_list)

    # compute only a subset of the divergence matrix
    bsc = bs(sts)
    div = np.zeros((num_indivs, num_targets))
    for k in range(num_targets):
        f = lambda x: [float(x[k]*(2-x[j])) for j in range(num_targets, num_targets + num_indivs)]
        div[:, k] = bsc.tree_stat_vector(simplified_node_list, weight_fun=f, windows=[0.0, ts.sequence_length])[0]

    div /= num_indivs
    div *= 1e-4

    dists = np.sqrt((locs[indivs, 0, np.newaxis] - locs[targets, 0].T) ** 2
                    + (locs[indivs, 1, np.newaxis] - locs[targets, 1].T) ** 2)

    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    fig = plt.figure(figsize=(6, 6 * ymax / xmax))
    colors = ['b', 'c', 'm', 'y']
    ax = fig.add_subplot(111)
    for k in range(num_targets):
        ax.scatter(dists[:,k], div[:,k],
                   s=10, 
                   c=colors[k])
    ax.set_xlabel("geographic distance")
    ax.set_ylabel("genetic distance")
    return fig


for script in ("flat_map.slim", "valleys.slim"):
    num_gens = 1
    treefile = sps.run_slim(script = script,
                            seed = 23, 
                            SIGMA = 0.4,
                            W = 50.0, 
                            NUMGENS = num_gens,
                            BURNIN = 10000)
    outbase = ".".join(treefile.split(".")[:-1])

    ts = pyslim.load(treefile)
    ts = sps.SpatialSlimTreeSequence(ts.recapitate(Ne = 1e3, recombination_rate = 1e-9), dim=2)

    locs = ts.individual_locations()
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    alive = ts.individuals_alive(0)
    left_third = np.logical_and(alive, locs[:, 0] < xmax / 3)
    right_third = np.logical_and(alive, locs[:, 0] > 2 * xmax / 3)
    targets = list(np.random.choice(np.where(left_third)[0], 2)) + \
                list(np.random.choice(np.where(right_third)[0], 2))

    fig = plot_isolation_by_distance(ts, targets, 100)
    plt.tight_layout() # please don't cut off my axis labels
    fig.savefig(outbase + ".ibd.pdf")

