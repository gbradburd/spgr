import pyslim, msprime
import numpy as np
import spatial_slim as sps
import os, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections as cs

usage = """
Usage:
    {} (script name)
""".format(sys.argv[0])

if len(sys.argv) != 2:
    raise ValueError(usage)

script = sys.argv[1]

def compute_ancestry(ts):
    """
    Find the mean proportion of the children's genome that each node is parental to
    """
    locs = ts.individual_locations()
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    alive = ts.individuals_alive(0)
    left_third = np.logical_and(alive, locs[:, 0] < xmax / 3)
    right_third = np.logical_and(alive, locs[:, 0] > 2 * xmax / 3)
    targets = list(np.random.choice(np.where(left_third)[0], 2)) + \
                list(np.random.choice(np.where(right_third)[0], 2))

    children_nodes = [ts.individual(ind).nodes for ind in targets]
    node_ancestry = ts.proportion_ancestry_nodes(children_nodes)
    return targets, node_ancestry

def plot_ancestry(ts, targets, times, node_ancestry):
    """
    Plot how the ancestry of an individual spreads out over space.
    Note that this doesn't look quite right because parents have the child's mass
    even *after* the child is born.
    """

    locs = ts.individual_locations()

    def size_fun(inds, k, scale=2000):
        return np.fromiter(map(lambda x: scale * sum(node_ancestry[k][ts.individual(x).nodes]), inds),
                           'float')

    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    figs = []
    for time in times:
        fig = plt.figure(figsize=(6, 6 * ymax / xmax))
        ax = fig.add_subplot(111)
        plt.axis('equal')
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, ymax)
        sps.plot_density(ts, time, ax, scatter=False)
        colors = sps.four_colors
        inds = ts.individuals_by_time(time)
        for k, child in enumerate(targets):
            ax.scatter(locs[inds, 0],
                       locs[inds, 1], 
                       sizes=size_fun(inds, k),
                       facecolor=colors[k],
                       edgecolor=None,
                       alpha=0.75)
        figs.append(fig)
    return figs

# for script in ("valleys.slim", "flat_map.slim"):
## pass in script on command-line

num_gens = 100
treefile = sps.run_slim(script = script,
                        seed = 23, 
                        SIGMA = 1.0,
                        W = 50.0, 
                        K = 5.0,
                        NUMGENS = num_gens)
outbase = ".".join(treefile.split(".")[:-1])

ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)


times = list(range(0, 100, 20)) + [ts.slim_generation - 1]

datafile = outbase + ".ancestry.txt"
if os.path.isfile(datafile):
    print(datafile, "already exists.")
    data = np.loadtxt(datafile)
    targets = np.array([np.int(u) for u in data[:, 0]])
    node_ancestry = data[:, 1:]
else:
    print(datafile, "does not exist, computing.")
    targets, node_ancestry = compute_ancestry(ts)
    data = np.column_stack([targets, node_ancestry])
    np.savetxt(datafile, data)

figs = plot_ancestry(ts, targets, times, node_ancestry)
for time, fig in zip(times, figs):
    fig.savefig(outbase + ".{:02d}.ancestry.pdf".format(time))
    plt.close(fig)
