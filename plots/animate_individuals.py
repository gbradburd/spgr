import pyslim, msprime
import numpy as np
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.collections as cs

def animate_individuals(ts, num_gens)
    # an animation of the individuals
    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(111)
    xmax = max([ind.location[0] for ind in ts.individuals()])
    ymax = max([ind.location[1] for ind in ts.individuals()])
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    # colors
    colormap = lambda x: plt.get_cmap("cool")(x/max(ts.individual_ages()))
    locs = ts.individual_locations()
    inds = ts.individuals_by_time(num_gens)
    next_inds = ts.individuals_by_time(num_gens - 1)
    circles = ax.scatter(locs[inds, 0], locs[inds, 1], s=10, 
                         edgecolors=colormap([0 for _ in inds]),
                         facecolors='none')
    filled = ax.scatter(locs[next_inds, 0], locs[next_inds, 1], s=10, 
                        facecolors=colormap([0 for _ in next_inds]),
                        edgecolors='none')
    lc = cs.LineCollection([], colors='black', linewidths=0.5)
    ax.add_collection(lc)

    def update(frame):
        inds = ts.individuals_by_time(frame)
        next_inds = ts.individuals_by_time(frame - 1)
        circles.set_offsets(locs[inds,:])
        filled.set_offsets(locs[next_inds,:])
        # color based on age so far
        circles.set_color(colormap(ts.individuals_age(frame)[inds]))
        filled.set_color(colormap(ts.individuals_age(frame)[next_inds]))
        if frame > 0:
            new_inds = inds[ts.individuals_age(frame)[inds] == 0]
            pcs = ts.get_individual_parents(new_inds, time=frame)
            lc.set_paths([locs[pc,:] for pc in pcs])
        return circles, filled, lc

    animation = ani.FuncAnimation(fig, update, 
                                  frames=np.linspace(num_gens, 1, num_gens))
    return animation

for script in ("flat_map.slim", "valleys.slim"):
    num_gens = 30
    treefile = sps.run_slim(script = script,
                            seed = 23, 
                            SIGMA = 0.25,
                            W = 50.0, 
                            NUMGENS = num_gens)
    outbase = ".".join(treefile.split(".")[:-1])

    ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

    animation = animate_individuals(ts, num_gens)
    animation.save(outbase + ".pop.mp4", writer='ffmpeg')

