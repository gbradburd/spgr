import pyslim, msprime
import numpy as np
import spatial_slim as sps
from msprime import BranchLengthStatCalculator as bs
import os, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

usage = """
Usage:
    {} (script name)
""".format(sys.argv[0])

if len(sys.argv) != 2:
    raise ValueError(usage)

script = sys.argv[1]

def draw_pie(ax, X, Y, ratios, size = 100, piecolors=["#ca0020", "#0571b0"]): 
    """
    Function to draw pie charts on map
    from https://stackoverflow.com/questions/51409257/exploding-wedges-of-pie-chart-when-plotting-them-on-a-map-python-matplotlib
    """
    xy = []
    s = []
    start = 0.0 
    for ratio in ratios:
        angles = np.linspace(2 * np.pi * start, 2 * np.pi * (start + ratio))
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()
        xy1 = np.column_stack([x, y])
        s1 = np.abs(xy1).max()
        xy.append(xy1)
        s.append(s1)
        start += ratio
    for k, (xyi, si) in enumerate(zip(xy,s)):
       ax.scatter([X], [Y], marker=xyi, s=size * si ** 2, edgecolor="k",  
                  facecolors=piecolors[k], linewidth=0.5, alpha=.75)


def compute_admixture(ts, time, admixture_time, num_targets):
    locs = ts.individual_locations()
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    alive = ts.individuals_alive(time)
    samples = np.random.choice(np.where(alive)[0], num_targets, replace=False)
    sample_nodes = ts.individual_nodes(samples, flatten=False)
    west_ancestors = np.logical_and(ts.individuals_alive(admixture_time),
                                    locs[:, 0] <= 1 * xmax / 2)
    east_ancestors = np.logical_and(ts.individuals_alive(admixture_time),
                                    locs[:, 0] > 1 * xmax / 2)
    west_ancestor_nodes = ts.individual_nodes(np.where(west_ancestors)[0])
    east_ancestor_nodes = ts.individual_nodes(np.where(east_ancestors)[0])

    def node_admixture(ts, sample_sets, ancestor_groups):
        """
        Find the proportion the genomes of each of sample_sets that inherits
        from each of ancestor_groups: if the output is D, then
        D[i,j] gives the proportion of sample_sets[i]'s ancestry that is contributed
        by ancestor_groups[j].
        """
        anc = ts.proportion_ancestry_nodes(sample_sets, show_progress=True)
        D = np.zeros((len(sample_sets), len(ancestor_groups)))
        for k, group in enumerate(ancestor_groups):
            D[:, k] = np.sum(anc[:, group], axis=1)
        return D

    admixture = node_admixture(ts, sample_nodes, [west_ancestor_nodes, east_ancestor_nodes])
    return samples, admixture


def plot_admixture(ts, time, samples, admixture):
    assert(len(samples) == admixture.shape[0])

    locs = ts.individual_locations()
    alive = ts.individuals_alive(time)
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])

    fig = plt.figure(figsize=(6, 6 * ymax / xmax))
    ax = fig.add_subplot(111)
    plt.axis('equal')
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    sps.plot_density(ts, time, ax)
    # ax.scatter(locs[alive, 0], locs[alive, 1],
    #            s=2, 
    #            alpha = 0.5,
    #            edgecolors='none',
    #            facecolors='black')
    for k, indiv in enumerate(samples):
        admix = admixture[k, :]
        admix /= sum(admix)
        draw_pie(ax, locs[indiv, 0], locs[indiv, 1], admix)

    return fig


num_gens = 301
treefile = sps.run_slim(script = script,
                        seed = 23, 
                        SIGMA = 4.0,
                        W = 50.0, 
                        K = 5.0,
                        NUMGENS = num_gens,
                        BURNIN=1)
outbase = ".".join(treefile.split(".")[:-1])
num_samples = 200
num_times = 5

# num_gens = 901
# treefile = sps.run_slim(script = script,
#                         seed = 23, 
#                         SIGMA = 3.0,
#                         W = 50.0, 
#                         K = 5.0,
#                         NUMGENS = num_gens,
#                         BURNIN=1)
# outbase = ".".join(treefile.split(".")[:-1])
# num_samples = 20
# num_times = 3

ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

for time in np.floor(np.linspace(0, num_gens - 1, num_times)):
    datafile = outbase + ".{}.admixture.txt".format(time)
    if os.path.isfile(datafile):
        print(datafile, "already exists.")
        data = np.loadtxt(datafile)
        samples = np.array([np.int(u) for u in data[:, 0]])
        admixture = data[:, 1:]
    else:
        print(datafile, "does not exist, computing.")
        samples, admixture = compute_admixture(ts, time, num_gens - 1, num_samples)
        data = np.column_stack([samples, admixture])
        np.savetxt(datafile, data)

    fig = plot_admixture(ts, time, samples, admixture)
    fig.savefig(outbase + ".{}.admixture.pdf".format(time))
    plt.close(fig)


