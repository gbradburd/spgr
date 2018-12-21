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

treefile = "flat_map.300.23.trees"
outbase = ".".join(treefile.split(".")[:-1])

ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

circle = ts.individuals_in_circle(center=[4,5], radius=2)
colors = ['red' if x else 'black' for x in np.isin(range(ts.num_individuals), circle)]

def plot_individuals(ts, outfile):
    # a snapshot of the individual locations
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.scatter(ts.individual_locations[:,0], ts.individual_locations[:,1],
               s=10, 
               c=colors)
    fig.savefig(outfile)


def animate_individuals(ts, outfile):
    # an animation of the individuals
    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(111)
    xmax = max([ind.location[0] for ind in ts.individuals()])
    ymax = max([ind.location[1] for ind in ts.individuals()])
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    # colors
    colormap = lambda x: plt.get_cmap("cool")(x/max(ts.individual_ages))
    locs = ts.individual_locations
    births = ts.individual_times
    inds = ts.individuals_by_time(ts.slim_generation)
    next_inds = ts.individuals_by_time(ts.slim_generation - 1)
    circles = ax.scatter(locs[inds, 0], locs[inds, 1], s=10, 
                         edgecolors=colormap([0 for _ in inds]),
                         facecolors='none')
    filled = ax.scatter(locs[next_inds, 0], locs[next_inds, 1], s=10, 
                        facecolors=colormap([0 for _ in inds]),
                        edgecolors='none')
    lc = cs.LineCollection([], colors='black', linewidths=0.5)
    ax.add_collection(lc)

    def update(frame):
        inds = ts.individuals_by_time(frame)
        next_inds = ts.individuals_by_time(frame - 1)
        circles.set_offsets(locs[inds,:])
        filled.set_offsets(locs[next_inds,:])
        # color based on age so far
        circles.set_color(colormap(births[inds] - frame))
        filled.set_color(colormap(births[next_inds] - frame))
        if frame > 0:
            new_inds = inds[births[inds] == frame]
            pcs = ts.get_individual_parents(new_inds, time=frame)
            lc.set_paths([locs[pc,:] for pc in pcs])
        return circles, filled, lc

    animation = ani.FuncAnimation(fig, update, 
                                  frames=np.linspace(ts.slim_generation, 1, ts.slim_generation))
    animation.save(outfile, writer='ffmpeg')


def animate_tree(ts, children, outfile):
    # an animation of the tree ancestral to an individual
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
    births = ts.individual_times
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
        circles.set_color(colormap(births[inds] - frame))
        newborns = children[births[children] == frame]
        pcs = ts.get_individual_parents(newborns)
        if len(pcs) > 0:
            children = np.concatenate((children, [a for a, b in pcs]))
            paths = paths + [locs[pc,:] for pc in pcs]
            lc.set_paths(paths)
        return circles, lc

    animation = ani.FuncAnimation(fig, update, 
                                  frames=np.linspace(0, ts.slim_generation, ts.slim_generation + 1))
    animation.save(outfile, writer='ffmpeg')

def animate_ancestry(ts, children, outfile):
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
    colormap = lambda x: plt.get_cmap("cool")(x/max(ts.individual_ages))
    locs = ts.individual_locations
    births = ts.individual_times
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
        circles.set_color(colormap(births[inds] - frame))
        circles.set_sizes(size_fun(inds)),
        return circles

    animation = ani.FuncAnimation(fig, update, 
                                  frames=np.linspace(0, ts.slim_generation - 1, ts.slim_generation))
    animation.save(outfile, writer='ffmpeg')


today = np.where(ts.individual_times == 0)[0]

if True:
    plot_individuals(ts, outbase + ".locations.pdf")
    animate_individuals(ts, outbase + ".pop.mp4")
    animate_tree(ts, np.random.choice(today, 10), outbase + ".trees.mp4")

animate_ancestry(ts, np.random.choice(today, 1), outbase + ".ancestry.mp4")
