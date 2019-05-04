import pyslim, msprime, tsinfer
import numpy as np
import spatial_slim as sps

import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.collections as cs


def animate_fitness(ts, num_gens):
    # an animation of the eventual genetic contributions of each individual.
    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(111)
    xmax = max([ind.location[0] for ind in ts.individuals()])
    ymax = max([ind.location[1] for ind in ts.individuals()])
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    # colors
    locs = ts.individual_locations
    # find the mean number of the present-day's genomes that each node is ancestral to
    today = ts.individuals_by_time(0)
    today_nodes = ts.individual_nodes(today, flatten=True)
    node_ancestry = len(today_nodes) * ts.proportion_ancestry_nodes([today_nodes])[0]
    indiv_ancestry = np.fromiter(map(lambda x: sum(node_ancestry[x.nodes]), ts.individuals()), 'float')

    def net_ancestry(time):
        # we need to subtract from each individual the values of any of their living children nodes: 
        # this is the only known situation when diploidy actually simplifies things!
        net_ancestry = indiv_ancestry.copy()
        indivs_alive = ts.individuals_alive(time)
        nodes_alive = indivs_alive[ts.tables.nodes.individual]
        nodes_alive[ts.tables.nodes.individual < 0] = False
        nodes = np.where(nodes_alive)[0]
        ins = [(ts.node(p).individual, c) for p,c in ts.node_children(nodes) if nodes_alive[c]]
        # testing
        ages = ts.individual_ages
        for i,n in ins:
            assert(ts.node(ts.individual(i).nodes[0]).time > time)
            assert(ts.node(ts.individual(i).nodes[0]).time - ages[i] <= time)
            assert(any(np.logical_and(np.isin(ts.tables.edges.parent, ts.individual(i).nodes),
                                      ts.tables.edges.child == n)))
        ch_dict = {n : set([ts.node(p).individual for p,c in ts.node_parents([n])]) for n in nodes}
        for c in ch_dict:
            assert(len(ch_dict[c]) == 1)
        ins_dict = {i : [c for p,c in ts.node_children(ts.individual(i).nodes)] for i in np.where(indivs_alive)[0]}
        for i in ins_dict:
            nn = ts.individual(i).nodes
            x = sum(np.isin(nn, today_nodes))
            for c in ins_dict[i]:
                x += node_ancestry[c]
            print(i, nn, x, indiv_ancestry[i], ins_dict[i])
            assert(x == indiv_ancestry[i])

        for ind, node in ins:
            if not (net_ancestry[ind] + 1e-8 > node_ancestry[node]):
                print(ind, node, net_ancestry[ind], node_ancestry[node])
            net_ancestry[ind] -= node_ancestry[node]
        return net_ancestry

    def size_fun(inds, scale=200):
        return np.fromiter(map(lambda x: 1 + scale * sum(node_ancestry[ts.individual(x).nodes]), inds),
                           'float')

    inds = ts.individuals_by_time(0)
    circles = ax.scatter(locs[inds, 0], locs[inds, 1], 
                         sizes=size_fun(inds),
                         color='black',
                         alpha=0.75)
    def update(frame):
        inds = ts.individuals_by_time(frame)
        circles.set_offsets(locs[inds,:])
        # color based on age so far
        circles.set_sizes(size_fun(inds))
        return circles

    start_gen = 0
    end_gen = num_gens - 1
    ngens = abs(end_gen - start_gen)
    frames = np.linspace(start_gen, end_gen, ngens)
    total_duration = 10 # seconds
    # interval is in milliseconds
    animation = ani.FuncAnimation(fig, update, 
                                  frames=frames,
                                  interval=total_duration * 1e3 / ngens)
    return animation


for script in ("flat_map.slim", "valleys.slim"):
    treefile = sps.run_slim(script = script,
                            seed = 23, 
                            SIGMA = 0.4,
                            W = 8.0, 
                            NUMGENS = 300)
    outbase = ".".join(treefile.split(".")[:-1])

    ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

    animation = animate_fitness(ts, num_gens)
    animation.save(outbase + ".fitness.mp4", writer='ffmpeg')

