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

def plot_ancestry(ts, children, times):
    """
    Plot how the ancestry of an individual spreads out over space.
    Note that this doesn't look quite right because parents have the child's mass
    even *after* the child is born.
    """

    locs = ts.individual_locations()
    # find the mean proportion of the children's genome that each node is parental to
    children_nodes = [ts.individual(ind).nodes for ind in children]
    node_ancestry = ts.proportion_ancestry_nodes(children_nodes)

    def size_fun(inds, k, scale=2000):
        return np.fromiter(map(lambda x: scale * sum(node_ancestry[k][ts.individual(x).nodes]), inds),
                           'float')

    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    figs = []
    for time in times:
        fig = plt.figure(figsize=(6, 6 * ymax / xmax))
        ax = fig.add_subplot(111)
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, ymax)
        colors = ['b', 'r', 'm', 'y']
        inds = ts.individuals_by_time(time)
        for k, child in enumerate(children):
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
today = np.where(ts.individual_times() == 0)[0]

locs = ts.individual_locations()
xmax = max(locs[:,0])
ymax = max(locs[:,1])
alive = ts.individuals_alive(0)
left_third = np.logical_and(alive, locs[:, 0] < xmax / 3)
right_third = np.logical_and(alive, locs[:, 0] > 2 * xmax / 3)
targets = list(np.random.choice(np.where(left_third)[0], 2)) + \
            list(np.random.choice(np.where(right_third)[0], 2))

times = list(range(0, 100, 20)) + [ts.slim_generation - 1]
figs = plot_ancestry(ts, targets, times)
for time, fig in zip(times, figs):
    fig.savefig(outbase + ".{:02d}.ancestry.pdf".format(time))
    plt.close(fig)
