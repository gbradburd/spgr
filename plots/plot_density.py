import pyslim, msprime
import numpy as np
import spatial_slim as sps
from msprime import BranchLengthStatCalculator as bs

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import scipy.stats

def plot_density(ts, time):
    """
    Plot a 2D kernel density estimate of the population density.
    """

    locs = ts.individual_locations()
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    alive = ts.individuals_alive(time)

    tlocs = locs[alive, :].T
    kde = scipy.stats.gaussian_kde(tlocs)
    X, Y = np.meshgrid(
            np.linspace(0.0, xmax, 51),
            np.linspace(0.0, ymax, 51))
    Z = kde([X.flatten(), Y.flatten()])
    Z.shape = X.shape

    fig = plt.figure(figsize=(6, 6 * ymax / xmax))
    ax = fig.add_subplot(111)
    ax.scatter(locs[alive, 0], locs[alive, 1],
               s=10,
               alpha=0.5,
               c='black',
               marker="o",
               edgecolors='none')
    ax.contour(X, Y, Z,
               colors='c',
               alpha=0.8)
    return fig


for script in ("valleys.slim", "flat_map.slim"):
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

    fig = plot_density(ts, 0)
    fig.savefig(outbase + ".density.pdf")
    plt.close(fig)


