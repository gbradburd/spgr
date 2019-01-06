import pyslim, msprime, tsinfer
import numpy as np
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.collections as cs

import importlib
importlib.reload(sps)


def animate_ancestry(ts, children, num_gens):
    # an animation of how the ancestry of an individual spreads out over space
    # note that this doesn't look quite right because parents have the child's mass
    # even *after* the child is born
    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(111)
    xmax = max([ind.location[0] for ind in ts.individuals()])
    ymax = max([ind.location[1] for ind in ts.individuals()])
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    # colors
    colormap = lambda x: plt.get_cmap("cool")(x/max(ts.individual_ages()))
    locs = ts.individual_locations()
    # find the mean proportion of the children's genome that each node is parental to
    child_nodes = [n for ind in children for n in ts.individual(ind).nodes]
    node_ancestry = ts.proportion_ancestry_nodes([child_nodes])[0]

    def size_fun(inds, scale=2000):
        return np.fromiter(map(lambda x: 1 + scale * sum(node_ancestry[ts.individual(x).nodes]), inds),
                           'float')

    inds = ts.individuals_by_time(0)
    circles = ax.scatter(locs[inds, 0], locs[inds, 1], 
                         sizes=size_fun(inds),
                         color=colormap([0 for _ in inds]),
                         alpha=0.75)
    def update(frame):
        inds = ts.individuals_by_time(frame)
        circles.set_offsets(locs[inds,:])
        # color based on age so far
        circles.set_color(colormap(ts.individuals_age(frame)[inds]))
        circles.set_sizes(size_fun(inds)),
        return circles

    animation = ani.FuncAnimation(fig, update, 
                                  frames=np.linspace(0, num_gens - 1, num_gens))
    return animation



for script in ("flat_map.slim", "valleys.slim"):
    num_gens = 300
    treefile = sps.run_slim(script = script,
                            seed = 23, 
                            SIGMA = 0.4,
                            W = 8.0, 
                            NUMGENS = num_gens)
    outbase = ".".join(treefile.split(".")[:-1])

    ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)
    today = np.where(ts.individual_times() == 0)[0]

    animation = animate_ancestry(ts, np.random.choice(today, 1), num_gens)
    animation.save(outbase + ".ancestry.mp4", writer='ffmpeg')
