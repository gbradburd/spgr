import pyslim, msprime
import numpy as np
import spatial_slim as sps
import os, sys

# we have the full pedigree for this long
num_gens = 200
treefile = sps.run_slim(script = "flat_map.slim",
                        seed = 23, 
                        SIGMA = 1.0,
                        W = 50.0, 
                        K = 5.0,
                        NUMGENS = num_gens,
                        BURNIN=1)
outbase = ".".join(treefile.split(".")[:-1])
# treefile = "valleys/run_1.0_50.0_5.0_100_23.trees"


ts = sps.SpatialSlimTreeSequence(pyslim.load(treefile), dim=2)

###
# Find effective generation time:
#   for each time t in the past,
#   take the average over all individuals alive at time t
#   of the average of their two parents' age at the time of their birth
#   where individuals are chosen proportionally to their genetic contribution to the modern generation
# So, as t increases, this should converge to the effective generation time.
#
# Said another way: 
#  1) pick a bit of today's genome
#  2) trace back to the individual, x, from whom it was inherited t time units in the past
#  3) pick a random one of x's two parents, whom we call y.
#  4) call A the age of y at the time x was born
#  5) We report the average value of A, as a function of t.

# find how much everyone contributed to today's nodes
alive = ts.individuals_alive(0)
alive_nodes = ts.individual_nodes(np.where(alive)[0])
node_ancestry = ts.proportion_ancestry_nodes([alive_nodes])[0]
ind_ancestry = np.zeros(ts.num_individuals)
for k, ind in enumerate(ts.individuals()):
    ind_ancestry[k] = sum(node_ancestry[ind.nodes])

# now find average parent age for everyone 
parents = ts.individual_parents_dict()
mean_parent_ages = np.zeros(ts.num_individuals) - 1.0
for ind in parents:
    ages = [ts.individual(x).time - ts.individual(ind).time for x in set(parents[ind])]
    if len(ages) == 2:
        mean_parent_ages[ind] = np.mean(ages)


# the max generation time is this long -
# we need to avoid looking at generations where we don't have their parents recorded
max_gen_time = 20
parent_age_by_time = np.zeros(num_gens - max_gen_time)
for t in range(num_gens - max_gen_time):
    alive = np.logical_and(ts.individuals_alive(t), mean_parent_ages > 0)
    probs = ind_ancestry * alive
    probs = probs / sum(probs)
    parent_age_by_time[t] = sum(probs * mean_parent_ages)


print("time_ago eff_generation_time")
for t in range(num_gens - max_gen_time):
    print("{} {}", t, parent_age_by_time[t])
