import pyslim, msprime
import numpy as np
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_individuals(ts, num_gens):
    # a snapshot of the individual locations
    xmax = max([ind.location[0] for ind in ts.individuals()])
    ymax = max([ind.location[1] for ind in ts.individuals()])
    fig = plt.figure(figsize=(9, 9 * ymax / xmax))
    ax = fig.add_subplot(111)
    plot_these = (ts.individual_times() <= num_gens)
    locs = ts.individual_locations()
    ax.scatter(locs[plot_these,0], locs[plot_these,1],
               s=10, 
               c='black')
    return fig

for script in ("flat_map.slim", "valleys.slim"):
    num_gens = 30
    treefile = sps.run_slim(script = script,
                            seed = 23, 
                            SIGMA = 0.25,
                            W = 50.0, 
                            NUMGENS = num_gens)
    outbase = ".".join(treefile.split(".")[:-1])

    ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

    fig = plot_individuals(ts, num_gens)
    fig.savefig(outbase + ".locations.pdf")
