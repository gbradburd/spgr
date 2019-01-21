import pyslim, msprime
import numpy as np
import spatial_slim as sps
from msprime import BranchLengthStatCalculator as bs
import os, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import scipy.stats

usage = """
Usage:
    {} (script name)
""".format(sys.argv[0])

if len(sys.argv) != 2:
    raise ValueError(usage)

script = sys.argv[1]

# for script in ("valleys.slim", "flat_map.slim"):
## pass in script on command-line

num_gens = 301
treefile = sps.run_slim(script = script,
                        seed = 23, 
                        SIGMA = 4.0,
                        W = 50.0, 
                        K = 5.0,
                        NUMGENS = num_gens,
                        BURNIN=1)
outbase = ".".join(treefile.split(".")[:-1])

ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)
locs = ts.individual_locations()
xmax = max(locs[:,0])
ymax = max(locs[:,1])

fig = plt.figure(figsize=(6, 6 * ymax / xmax))
ax = fig.add_subplot(111)
plt.axis('equal')
ax.set_xlim(0, xmax)
ax.set_ylim(0, ymax)
sps.plot_density(ts, 0, ax)
fig.savefig(outbase + ".density.pdf")
plt.close(fig)


