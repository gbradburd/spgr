import pyslim, msprime
import numpy as np
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.collections as cs
import matplotlib.patches as ptch

def plot_flux(ts, inout_fns, plot_ngens):
    """
    Plot arrows between all parent-offspring pairs for which each fn in
    inout_fns applied to their two locations are either True, False (in: red
    arrow) or False, True (out: blue arrow). The functions in the list
    `inout_fns` should take an array of locations and return a logical vector.
    """
    xmax = max([ind.location[0] for ind in ts.individuals()])
    ymax = max([ind.location[1] for ind in ts.individuals()])
    fig = plt.figure(figsize=(9, 9 * ymax / xmax))
    ax = fig.add_subplot(111)
    plt.axis('equal')
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)

    locs = ts.individual_locations()
    births = ts.individual_times()
    deaths = births - ts.individual_ages()
    alive = (deaths <= plot_ngens)  # alive during the last ngens 
    ips = ts.individual_parents(list(range(ts.num_individuals)))
    born = (births[ips[:, 1]] <= plot_ngens) # pairs born during the last ngens
    # set up plot with individual locations
    circles = ax.scatter(locs[alive, 0], locs[alive, 1], s=10, 
                         edgecolors='none', 
                         alpha=0.5,
                         facecolors='c')

    for inout_fn in inout_fns:
        parents_in = inout_fn(locs[ips[:, 0], :])
        children_in = inout_fn(locs[ips[:, 1], :])
        in_pairs = np.logical_and(born,
                        np.logical_and(parents_in, 
                                       np.logical_not(children_in)))
        out_pairs = np.logical_and(born,
                        np.logical_and(np.logical_not(parents_in), 
                                       children_in))
        assert(sum(np.logical_and(in_pairs, out_pairs)) == 0)

        for pairs, color in ((in_pairs, 'r'), (out_pairs, 'b')):
            if sum(pairs) > 0:
                ax.quiver(locs[ips[pairs, 0], 0],  # X
                          locs[ips[pairs, 0], 1],  # X
                          locs[ips[pairs, 1], 0]
                           - locs[ips[pairs, 0], 0], # dX
                          locs[ips[pairs, 1], 1]
                           - locs[ips[pairs, 0], 1], # dY
                          color = color,
                          alpha = 0.25,
                          width = 0.03,
                          units='xy', scale=1)
    return fig

treefiles = [
    sps.run_slim(script = "valleys.slim",
                 seed = 23, 
                 SIGMA = 1.5,
                 W = 50.0, 
                 K = 5.0,
                 NUMGENS = 30),
    sps.run_slim(script = "flat_map.slim",
                 seed = 23, 
                 SIGMA = 0.4,
                 W = 50.0, 
                 K = 5.0,
                 NUMGENS = 30)
    ]

for treefile in treefiles:
    outbase = ".".join(treefile.split(".")[:-1])

    ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

    # number of time steps to plot the flux for
    plot_ngens = 10

    xmax = max([ind.location[0] for ind in ts.individuals()])
    ymax = max([ind.location[1] for ind in ts.individuals()])
    line_x = [xmax / 3., xmax * 2. / 3.]

    def make_line_fn(x):
        def f(locs):
            return (locs[:, 0] > x)
        return f

    line_inout_fns = [make_line_fn(x) for x in line_x]

    fig = plot_flux(ts, line_inout_fns, plot_ngens)
    ax = fig.axes[0]
    for x in line_x:
        ax.axvline(x, color='black', alpha=0.5, dashes=[3, 1], linewidth=2.0)
    fig.savefig(outbase + ".line_flux.pdf")

    circle_xy = (xmax / 2, ymax / 2)
    circle_rad = xmax / 3

    def circle_inout(locs):
        return (np.square(locs[:, 0] - circle_xy[0])
                 + np.square(locs[:, 1] - circle_xy[1]) < circle_rad**2)

    fig = plot_flux(ts, [circle_inout], plot_ngens)
    ax = fig.axes[0]
    circle = ptch.Circle(circle_xy, circle_rad,
                         edgecolor='black', alpha=0.5, linestyle=(0., (3., 1)), linewidth=2.0)
    ax.add_patch(circle)
    fig.savefig(outbase + ".circle_flux.pdf")
