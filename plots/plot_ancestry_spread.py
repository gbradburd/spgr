import pyslim, msprime
import numpy as np
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections as cs


def plot_ancestry(ts, children, time):
    """
    Plot how the ancestry of an individual spreads out over space.
    Note that this doesn't look quite right because parents have the child's mass
    even *after* the child is born.
    """

    inds = ts.individuals_by_time(time)
    locs = ts.individual_locations()
    # find the mean proportion of the children's genome that each node is parental to
    child_nodes = [n for ind in children for n in ts.individual(ind).nodes]
    node_ancestry = ts.proportion_ancestry_nodes([child_nodes])[0]

    def size_fun(inds, scale=2000):
        return np.fromiter(map(lambda x: 1 + scale * sum(node_ancestry[ts.individual(x).nodes]), inds),
                           'float')

    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    fig = plt.figure(figsize=(6, 6 * ymax / xmax))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    circles = ax.scatter(locs[inds, 0], locs[inds, 1], 
                         sizes=size_fun(inds),
                         color='c',
                         alpha=0.75)
    circles.set_offsets(locs[inds,:])
    circles.set_sizes(size_fun(inds))

    return fig


for script in ("flat_map.slim", "valleys.slim"):
    num_gens = 300
    treefile = sps.run_slim(script = script,
                            seed = 23, 
                            SIGMA = 0.4,
                            W = 50.0, 
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

    for time in range(0, 100, 5):
        fig = plot_ancestry(ts, targets, time)
        fig.savefig(outbase + ".{:02d}.ancestry.pdf".format(time))
        plt.close(fig)
