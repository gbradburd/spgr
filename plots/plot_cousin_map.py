import pyslim, msprime
import numpy as np
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.collections as cs
import matplotlib.patches as ptch

import scipy.sparse as sparse

np.random.seed(32)

def plot_cousins(ts, focal_individuals, max_hops = 9):
    """
    For each extant individual, plot circle sizes proportional to
    the closest degree of relatedness they have to a focal individual,
    colored by individual.
    """
    xmax = max([ind.location[0] for ind in ts.individuals()])
    ymax = max([ind.location[1] for ind in ts.individuals()])
    fig = plt.figure(figsize=(9, 9 * ymax / xmax))
    ax = fig.add_subplot(111)
    plt.axis('equal')
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)

    relnesses = [ts.relatedness(ts.individual(u).nodes, max_hops) for u in focal_individuals]
    locs = ts.individual_locations()
    alive = ts.individuals_alive(0)
    # base locations
    circles = ax.scatter(locs[alive, 0], locs[alive, 1], 
                         s=2,
                         edgecolors='none', 
                         alpha=0.5,
                         facecolors='c')
    colors = ['g', 'r', 'm', 'y']
    for x, col in zip(relnesses, colors):
        ind_x = np.repeat(np.inf, ts.num_individuals)
        for u in np.where(x)[0]:
            ind_x[ts.node(u).individual] = min(ind_x[ts.node(u).individual], x[u])
        ind_x[ind_x == np.inf] = np.nan
        circles = ax.scatter(locs[alive, 0], locs[alive, 1], 
                             s=1000 * (2 ** (- ind_x[alive])), 
                             edgecolors='none', 
                             alpha=0.5,
                             facecolors=col)

    return fig


treefiles = [
    sps.run_slim(script = "valleys.slim",
                 seed = 23, 
                 SIGMA = 1.0,
                 W = 50.0, 
                 K = 5.0,
                 NUMGENS = 100),
    sps.run_slim(script = "flat_map.slim",
                 seed = 23, 
                 SIGMA = 1.0,
                 W = 50.0, 
                 K = 5.0,
                 NUMGENS = 100)
    ]

for treefile in treefiles:
    outbase = ".".join(treefile.split(".")[:-1])

    ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

    locs = ts.individual_locations()
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    alive = ts.individuals_alive(0)
    left_third = np.logical_and(alive, locs[:, 0] < xmax / 3)
    right_third = np.logical_and(alive, locs[:, 0] > 2 * xmax / 3)
    focal_individuals = list(np.random.choice(np.where(left_third)[0], 2)) + \
                         list(np.random.choice(np.where(right_third)[0], 2))

    fig = plot_cousins(ts, focal_individuals)
    fig.savefig(outbase + ".cousins.pdf")
