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

def compute_isolation_by_distance(ts, targets, num_indivs):
    """
    Compute mean sequence divergence of num_individuals to each of the targets
    """
    # a snapshot of the individual locations
    locs = ts.individual_locations()
    alive = ts.individuals_alive(0)
    other_alive_inds = np.array(list(set(np.where(alive)[0]) - set(targets)))
    indivs = list(np.random.choice(other_alive_inds, num_indivs, replace=False))
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

    dists = np.sqrt(np.power(locs[indivs, 0, np.newaxis] - locs[targets, 0].T, 2)
                    + np.power(locs[indivs, 1, np.newaxis] - locs[targets, 1].T, 2))

    return div, dists


def plot_isolation_by_distance(ts, div, dists):
    num_targets = div.shape[1]
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    fig = plt.figure(figsize=(6, 6 * ymax / xmax))
    ax = fig.add_subplot(111)
    if True:
        # in order, with different markers
        for k in range(num_targets):
            ax.scatter(dists[:, k], div[:, k],
                       s = 20, 
                       alpha = 0.75,
                       edgecolors = 'none',
                       facecolors = sps.four_colors[k],
                       marker=sps.four_markers[k])
    else:
        # randomized order
        colors = np.repeat(sps.four_colors, div.shape[0])
        x = dists[:, 0]
        y = div[:, 0]
        for k in range(1, num_targets):
            x = np.concatenate([x, dists[:, k]])
            y = np.concatenate([y, div[:, k]])
        po = np.random.choice(len(x), len(x), replace=False)
        ax.scatter(x[po], y[po],
                   s = 10, 
                   alpha = 0.5,
                   marker = ".",
                   c = colors[po])
    ax.set_xlabel("geographic distance")
    ax.set_ylabel("genetic distance")
    return fig


# for script in ("valleys.slim", "flat_map.slim"):
## pass in script on command-line

num_gens = 1
treefile = sps.run_slim(script = script,
                        seed = 23, 
                        SIGMA = 0.4,
                        W = 50.0, 
                        NUMGENS = num_gens,
                        BURNIN = 10000)
outbase = ".".join(treefile.split(".")[:-1])

recapfile = outbase + ".recap.trees"
if os.path.isfile(recapfile):
    ts = sps.SpatialSlimTreeSequence(pyslim.load(recapfile), dim=2)
else:
    ts = pyslim.load(treefile)
    ts = sps.SpatialSlimTreeSequence(ts.recapitate(Ne = 1e3, recombination_rate = 1e-9), dim=2)
    ts.dump(recapfile)

num_targets = 4
num_indivs = 200
locs = ts.individual_locations()
xmax = max(locs[:,0])
ymax = max(locs[:,1])
alive = ts.individuals_alive(0)
left_third = np.logical_and(alive, locs[:, 0] < xmax / 3)
right_third = np.logical_and(alive, locs[:, 0] > 2 * xmax / 3)
targets = list(np.random.choice(np.where(left_third)[0], int(num_targets / 2))) + \
            list(np.random.choice(np.where(right_third)[0], int(num_targets / 2)))

datafile = outbase + ".ibd.txt"
if os.path.isfile(datafile):
    print(datafile, "already exists.")
    data = np.loadtxt(datafile)
    div = data[:, :num_targets]
    dists = data[:, num_targets:]
else:
    print(datafile, "does not exist, computing.")
    div, dists = compute_isolation_by_distance(ts, targets, num_indivs)
    data = np.column_stack([div, dists])
    np.savetxt(datafile, data)

fig = plot_isolation_by_distance(ts, div, dists)
plt.tight_layout() # please don't cut off my axis labels
fig.savefig(outbase + ".ibd.pdf")

