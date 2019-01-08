import pyslim, msprime
import numpy as np
import spatial_slim as sps
from msprime import BranchLengthStatCalculator as bs

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_isolation_by_distance(ts, num_indivs):
    # a snapshot of the individual locations
    locs = ts.individual_locations()
    alive = ts.individuals_alive(0)
    indivs = np.random.choice(np.where(alive)[0], num_indivs, replace=False)

    # simplify so this goes faster
    individual_node_list = []
    simplified_node_list = []
    k = 0
    for ind in indivs:
        individual_node_list.extend(ts.individual(ind).nodes)
        simplified_node_list.append([k, k+1])
        k = k + 2

    sts = ts.simplify(individual_node_list)
    bsc = bs(sts)
    div = bsc.divergence_matrix(simplified_node_list, [0.0, ts.sequence_length])

    dists = np.sqrt((locs[indivs, 0, np.newaxis] - locs[indivs, 0].T) ** 2
                    + (locs[indivs, 1, np.newaxis] - locs[indivs, 1].T) ** 2)
    ut = np.triu_indices(len(indivs), k=1)

    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    fig = plt.figure(figsize=(6, 6 * ymax / xmax))
    ax = fig.add_subplot(111)
    ax.scatter(dists[ut], div[0][ut],
               s=30, 
               c='black')
    return fig


def plot_heterozygosity(ts):
    # DO NOT DO THIS: subset and simplify first
    het = [ts.pairwise_diversity(ts.individual(k).nodes) 
            for k in np.where(ts.individuals_alive(0))[0]]



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

    fig = plot_isolation_by_distance(ts, 200)
    fig.savefig(outbase + ".ibd.pdf")

