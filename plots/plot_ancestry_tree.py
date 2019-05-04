import pyslim, msprime
import numpy as np
import spatial_slim as sps
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.collections as cs

usage = """
Usage:
    {} (script name)
""".format(sys.argv[0])

if len(sys.argv) != 2:
    raise ValueError(usage)

script = sys.argv[1]


def animate_tree(ts, children, num_gens):
    """
    An animation of the tree ancestral to an individual.
    """
    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(111)
    xmax = max([ind.location[0] for ind in ts.individuals()])
    ymax = max([ind.location[1] for ind in ts.individuals()])
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    # colors
    colormap = lambda x: plt.get_cmap("cool")(x/max(ts.individual_ages))
    # treecolors = [plt.get_cmap("viridis")(x) for x in np.linspace(0, 1, len(children))]
    locs = ts.individual_locations
    inds = ts.individuals_by_time(0)
    circles = ax.scatter(locs[inds, 0], locs[inds, 1], s=10, 
                         edgecolors=colormap([0 for _ in inds]),
                         facecolors='none')
    paths = []
    lc = cs.LineCollection(paths, linewidths=0.5)
    ax.add_collection(lc)

    def update(frame):
        nonlocal children
        nonlocal paths
        inds = ts.individuals_by_time(frame)
        circles.set_offsets(locs[inds,:])
        # color based on age so far
        circles.set_color(colormap(ts.individuals_age(frame)[inds]))
        newborns = children[ts.individuals_age(frame)[children] == 0]
        pcs = ts.individual_parents(newborns)
        if len(pcs) > 0:
            children = np.concatenate((children, pcs[:,0]))
            paths = paths + [locs[pc,:] for pc in pcs]
            lc.set_paths(paths)
        return circles, lc

    animation = ani.FuncAnimation(fig, update, 
                                  frames=np.linspace(0, num_gens, num_gens + 1))
    return animation


num_gens = 300
treefile = sps.run_slim(script = script,
                        seed = 23, 
                        SIGMA = 0.4,
                        W = 8.0, 
                        NUMGENS = num_gens)
outbase = ".".join(treefile.split(".")[:-1])

ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

today = np.where(ts.individual_times == 0)[0]
animation = animate_tree(ts, np.random.choice(today, 10), num_gens)
animation.save(outbase + ".trees.mp4", writer='ffmpeg')

