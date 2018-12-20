import pyslim, msprime
import numpy as np


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

    def individuals_in_circle(self, center, radius):
        """
        Returns the IDs of individuals at distance less than or equal to `radius` of `center.

        :param float center: The coordinates of the center.
        :param float radius: The radius of the circle.
        """
        dists = self.individual_distance_to_point(center)
        return np.where(dists <= radius)[0]

    def individuals_alive(self, time):
        """
        Returns a logical array of length equal to the number of individuals that
        indicates whether the given individual was alive at the given time.

        :param float time: The time ago.
        """
        births = self.individual_times
        ages = self.individual_ages
        return np.logical_and(births >= time, births - ages <= time)

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
        if time is not None:
            alive = self.individuals_alive(time)
        else:
            alive = np.repeat(True, self.num_individuals)
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
