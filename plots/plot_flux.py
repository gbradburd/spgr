import pyslim, msprime
import numpy as np
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.collections as cs

import importlib
importlib.reload(sps)

treefile = "flat_map.300.23.trees"
outbase = ".".join(treefile.split(".")[:-1])

ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

def plot_flux(ts, inout_fn):
    """
    Plot arrows between all parent-offspring pairs for which inout_fn() applied
    to their two locations are either True, False (in: red arrow) or False, True
    (out: blue arrow). The function `inout_fn` should take an array of locations
    and return a logical vector.
    """
    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(111)
    plt.axis('equal')
    xmax = max([ind.location[0] for ind in ts.individuals()])
    ymax = max([ind.location[1] for ind in ts.individuals()])
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)

    locs = ts.individual_locations()
    ips = ts.get_individual_parents(list(range(ts.num_individuals)))
    parents_in = inout_fn(locs[ips[:, 0], :])
    children_in = inout_fn(locs[ips[:, 1], :])
    in_pairs = np.logical_and(parents_in, np.logical_not(children_in))
    out_pairs = np.logical_and(np.logical_not(parents_in), children_in)

    circles = ax.scatter(locs[:, 0], locs[:, 1], s=10, 
                         edgecolors='none',
                         facecolors='c')
    for pairs, color in ((in_pairs, 'r'), (out_pairs, 'b')):
        if sum(pairs) > 0:
            ax.quiver(locs[ips[pairs, 0], 0],  # X
                      locs[ips[pairs, 0], 1],  # X
                      locs[ips[pairs, 1], 0]
                       - locs[ips[pairs, 0], 0], # dX
                      locs[ips[pairs, 1], 1]
                       - locs[ips[pairs, 0], 1], # dY
                      color = color,
                      units='xy', scale=1)

    return fig


xmax = max([ind.location[0] for ind in ts.individuals()])
ymax = max([ind.location[1] for ind in ts.individuals()])
line_x = xmax / 2

def line_inout(locs):
    return locs[:, 0] > line_x


fig = plot_flux(ts, line_inout) 
fig.savefig(outbase + ".line_flux.png")

circle_xy = (xmax / 2, ymax / 2)
circle_rad = xmax / 5

def circle_inout(locs):
    return (np.square(locs[:, 0] - circle_xy[0])
             + np.square(locs[:, 1] - circle_xy[1]) < circle_rad)

fig = plot_flux(ts, circle_inout)
fig.savefig(outbase + ".circle_flux.png")
