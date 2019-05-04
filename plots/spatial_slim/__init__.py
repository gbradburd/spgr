import os, subprocess
import scipy.stats
import numpy as np
import matplotlib

from .spatial_slim_ts import *

# pop_width = 8.0, 
# numgens = 100, 
# sigma = 1.0, 
# K = 5.0,

def run_slim(script, seed = 23, 
             **kwargs):
    scriptbase = "_".join(script.split(".")[:-1])
    if not os.path.isdir(scriptbase):
        os.mkdir(scriptbase)
    base = os.path.join(scriptbase, 
            "_".join(["run"] + list(map(str, kwargs.values())) + [str(seed)]))
    treefile = base + ".trees"
    if os.path.isfile(treefile):
        print(treefile, "already exists.")
    else:
        logfile = base + ".log"
        slim_command = ["slim", "-s {}".format(seed)]
        slim_command += ["-d {}={}".format(k, v) for k, v in kwargs.items()]
        slim_command += ["-d \"OUTPATH='{}'\"".format(treefile), script]
        print(" ".join(slim_command))
        with open(logfile, "w") as log:
            subprocess.call(" ".join(slim_command), shell=True, stdout=log)
        if not os.path.isfile(treefile):
            raise ValueError("SLiM failed to produce output file {}".format(treefile))
    return(treefile)

# four_colors = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"]
four_colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a"]  # Dark2
four_markers = ["v", "<", "^", ">"]


def plot_density(ts, time, ax, scatter=True, alpha=0.8):
    """
    Plot a 2D kernel density estimate of the population density.
    """

    locs = ts.individual_locations[:,:2]
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    alive = ts.individuals_alive(time)
    tlocs = locs.T
    kde = scipy.stats.gaussian_kde(tlocs)
    X, Y = np.meshgrid(
            np.linspace(0.0, xmax, 51),
            np.linspace(0.0, ymax, 51))
    Z = kde([X.flatten(), Y.flatten()])
    Z.shape = X.shape
    if scatter:
        ax.scatter(locs[alive, 0], locs[alive, 1],
                   s=10,
                   alpha=0.5,
                   c='black',
                   marker="o",
                   edgecolors='none')
    ax.contour(X, Y, Z,
               colors='c',
               alpha=alpha,
               zorder=-1)

