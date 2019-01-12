import pyslim, msprime
import numpy as np
import spatial_slim as sps
from msprime import BranchLengthStatCalculator as bs

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def draw_pie(ax, X, Y, ratios, size = 200, piecolors=['c', 'm']): 
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


def plot_admixture(ts, time, admixture_time, num_targets):

    locs = ts.individual_locations()
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    alive = ts.individuals_alive(time)
    samples = np.random.choice(np.where(alive)[0], num_targets, replace=False)
    sample_nodes = ts.individual_nodes(np.where(samples)[0], flatten=False)
    west_ancestors = np.logical_and(ts.individuals_alive(admixture_time),
                                    locs[:, 0] <= 2 * xmax / 3)
    east_ancestors = np.logical_and(ts.individuals_alive(admixture_time),
                                    locs[:, 0] > 2 * xmax / 3)
    west_ancestor_nodes = ts.individual_nodes(np.where(west_ancestors)[0])
    east_ancestor_nodes = ts.individual_nodes(np.where(east_ancestors)[0])

    def node_admixture(ts, nodes, sample_sets):
        """
        Find the amount of each of nodes' genome that inherits from each of sample_sets.
        """
        anc = ts.proportion_ancestry_nodes(nodes, show_progress=True)
        N = len(nodes)
        K = len(sample_sets)
        D = np.zeros((N, K))
        for k, sample_set in enumerate(sample_sets):
            D[:, k] = np.sum(anc[:, sample_set], axis=1)

        return D

    admixture = node_admixture(ts, sample_nodes, [west_ancestor_nodes, east_ancestor_nodes])

    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    fig = plt.figure(figsize=(6, 6 * ymax / xmax))
    colors = ['b', 'c', 'm', 'y']
    ax = fig.add_subplot(111)
    ax.scatter(locs[alive, 0], locs[alive, 1],
               s=10, 
               c='black')
    for k, indiv in enumerate(samples):
        admix = admixture[k, :]
        admix /= sum(admix)
        print(locs[indiv, :], k, admix)
        draw_pie(ax, locs[indiv, 0], locs[indiv, 1], admix)

    return fig


for script in ("valleys.slim", "flat_map.slim"):
    num_gens = 100
    treefile = sps.run_slim(script = script,
                            seed = 23, 
                            SIGMA = 1.0,
                            W = 50.0, 
                            K = 5.0,
                            NUMGENS = num_gens)
    outbase = ".".join(treefile.split(".")[:-1])

    ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

    for time in [0, 25, 50, 75, 95]:
        fig = plot_admixture(ts, time, 99, 100)
        fig.savefig(outbase + ".{}.admixture.pdf".format(time))
        plt.close(fig)


