import pyslim, msprime
import numpy as np
import tqdm


class SpatialSlimTreeSequence(pyslim.SlimTreeSequence):

    def __init__(self, ts, dim=3):
        super(SpatialSlimTreeSequence, self).__init__(ts)
        self.dim = dim

    def individual_distance_to_point(self, point):
        """
        Returns the array of distances of each individual's location to the point
        whose coordinates are given by `point`.

        :param float point: The coordinates of the point. If this has length less than
            three, remaining coordinates will be set to 0.0.
        """
        point = np.concatenate((np.array(point), np.zeros(self.dim - len(point))))
        if point.shape != (self.dim,):
            raise ValueError("point not of the correct shape: must be coercible to"
                             + "a vector of length self.dim or less.")

        return np.sqrt(np.sum((self.individual_locations - point) ** 2, axis=1))

    def individuals_in_circle(self, center, radius, time=None):
        """
        Returns the IDs of individuals at distance less than or equal to `radius` of `center.

        :param float center: The coordinates of the center.
        :param float radius: The radius of the circle.
        """
        dists = self.individual_distance_to_point(center)
        these = np.logical_and(dists <= radius, self.individuals_alive(time))
        return np.where(these)[0]

    def individuals_alive(self, time):
        """
        Returns a logical array of length equal to the number of individuals that
        indicates whether the given individual was alive at the given time.

        :param float time: The time ago. If `time` is None, then everyone will be
            counted as being alive.
        """
        if time is None:
            out = np.repeat(True, self.num_individuals)
        else:
            births = self.individual_times
            ages = self.individual_ages
            out = np.logical_and(births >= time, births - ages <= time)
        return out

    def individuals_by_time(self, time):
        """
        Returns the IDs of individuals alive at the given time ago. Note that `age` is
        recorded while the individual is still alive, and starts at 0. TODO: clarify
        whether this is at `early()` or `late()` or what.

        :param float time: The amount of time ago.
        """
        return np.where(self.individuals_alive(time))[0]

    def get_node_parents(self, children, left=0.0, right=None):
        """
        Returns a list of all (parent, child) pairs of node IDs for the given 
        children, such that child inherited from parent somewhere in the region
        [left, right).

        :param int children: The node IDs of the children.
        :param float left: The left end of the portion of genome considered.
        :param float right: The right end of the portion of genome considered.
            Defaults to the sequence length.
        """
        if len(children) == 0:
            return []
        if max(children) >= self.num_nodes or min(children) < 0:
            raise ValueError("Node child index out of bounds.")
        if right is None:
            right = self.sequence_length
        if left < 0 or right > self.sequence_length or left > right:
            raise ValueError("Illegal left, right bounds.")
        edges = self.tables.edges
        yesthese = np.logical_and(np.isin(edges.child, children),
                                          edges.left < right,
                                          edges.right >= left)
        return zip(edges.parent[yesthese], edges.child[yesthese])

    def get_individual_parents(self, children, time=None, left=0.0, right=None):
        """
        Returns a list of all (parent, child) pairs of individual IDs for the given 
        children, such that child inherited from parent somewhere in the region
        [left, right), and the parent is alive at the given time.

        :param int children: The individual IDs of the children.
        :param float time: The time ago the parent should be alive. Defaults to
            no constraint.
        :param float left: The left end of the portion of genome considered.
            Defaults to zero.
        :param float right: The right end of the portion of genome considered.
            Defaults to the sequence length.
        """
        if len(children) == 0:
            return []
        if max(children) >= self.num_individuals or min(children) < 0:
            raise ValueError("Individual child index out of bounds.")
        child_nodes = [x for y in children for x in self.individual(y).nodes]
        node_parents = self.get_node_parents(child_nodes, left=left, right=right)
        alive = self.individuals_alive(time)
        out = [(self.node(a).individual, self.node(b).individual) for a, b in node_parents]
        return [(a, b) for (a, b) in out 
                if a is not msprime.NULL_INDIVIDUAL and b is not msprime.NULL_INDIVIDUAL
                   and alive[a]]
    @property
    def individual_times(self):
        return np.fromiter(map(lambda x: x.time, self.individuals()), 'float')

    @property
    def individual_ages(self):
        return np.fromiter(map(lambda x: pyslim.decode_individual(x.metadata).age, self.individuals()), 'int')

    @property
    def individual_populations(self):
        return np.fromiter(map(lambda x: pyslim.decode_individual(x.metadata).population, self.individuals()), 'int')

    @property
    def individual_locations(self):
        locations = self.tables.individuals.location
        locations.shape = (int(len(locations)/3), 3)
        return locations[:,:self.dim]

    def proportion_ancestry_nodes(self, sample_sets, show_progress=False):
        """
        Computes for each node the proportion of the genomes in each of sample sets
        inheriting from that node.
        """
        # Check the inputs (could be done more efficiently here)
        all_samples = set()
        for sample_set in sample_sets:
            U = set(sample_set)
            if len(U) != len(sample_set):
                raise ValueError("Cannot have duplicate values within set")
            if len(all_samples & U) != 0:
                # TODO: is this necessary? (shouldn't be)
                raise ValueError("Sample sets must be disjoint")
            all_samples |= U

        K = len(sample_sets)
        A = np.zeros((K, self.num_nodes))
        parent = np.zeros(self.num_nodes, dtype=int) - 1
        sample_count = np.zeros((K, self.num_nodes), dtype=int)
        last_update = np.zeros(self.num_nodes)
        total_length = np.zeros(self.num_nodes)

        def update_counts(edge, sign):
            # Update the counts and statistics for a given node. Before we change the
            # node counts in the given direction, check to see if we need to update
            # statistics for that node. When a node count changes, we add the
            # accumulated statistic value for the span since that node was last updated.
            v = edge.parent
            while v != -1:
                if last_update[v] != left:
                    total = np.sum(sample_count[:, v])
                    if total != 0:
                        length = left - last_update[v]
                        for j in range(K):
                            A[j, v] += length * sample_count[j, v]
                    last_update[v] = left
                for j in range(K):
                    sample_count[j, v] += sign * sample_count[j, edge.child]
                v = parent[v]

        # Set the intitial conditions.
        for j in range(K):
            for u in sample_sets[j]:
                sample_count[j][u] = 1

        progress_iter = tqdm.tqdm(
            self.edge_diffs(), total=self.num_trees, disable=not show_progress)
        for (left, right), edges_out, edges_in in progress_iter:
            for edge in edges_out:
                parent[edge.child] = -1
                update_counts(edge, -1)
            for edge in edges_in:
                parent[edge.child] = edge.parent
                update_counts(edge, +1)

        # Finally, add the stats for the last tree and normalise by the total
        # length that each node was an ancestor to > 0 samples.
        for v in range(self.num_nodes):
            total = np.sum(sample_count[:, v])
            if total != 0:
                length = self.sequence_length - last_update[v]
                for j in range(K):
                    A[j, v] += length * sample_count[j, v]
        # renormalize
        for j in range(K):
            A[j,:] /= len(sample_sets[j]) * self.sequence_length
        return A

    def proportion_ancestry_nodes_slow(self, sample_set, left=0.0, right=None):
        """
        Computes an array of values that gives, for each node, the total amount of
        genome in `sample_set` that descent from that node.

        IS VERY SLOOOOOW

        :param int sample_set: The children.
        :param float time: The time ago the parent should be alive. Defaults to
            no constraint.
        :param float left: The left end of the portion of genome considered.
            Defaults to zero.
        :param float right: The right end of the portion of genome considered.
            Defaults to the sequence length.
        """
        n = len(sample_set)
        if right is None:
            right = self.sequence_length
        if left < 0 or right > self.sequence_length or left > right:
            raise ValueError("Illegal left, right bounds.")
        trees = self.trees(sample_counts=True, tracked_leaves=sample_set)
        # first count up how much shared ancestry each node has
        ancestry = np.zeros(self.num_nodes)
        for tree in trees:
            if tree.interval[0] > right:
                break
            while tree.interval[1] < left:
                next
            tree_span = min(right, tree.interval[1]) - max(left, tree.interval[0])
            for n in tree.nodes():
                ancestry[n] += tree.num_tracked_samples(n) * tree_span
        ancestry /= len(sample_set) * ts.sequence_length
        return ancestry

