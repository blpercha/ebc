from collections import defaultdict
import random
import sys
import math

import numpy as np
from numpy.ma import divide, outer, sqrt
from numpy.random.mtrand import random_sample

from matrix import SparseMatrix

INFINITE = 1e10


class EBC:
    def __init__(self, matrix, n_clusters, max_iterations=10, jitter_max=1e-10, objective_tolerance=0.01):
        """ To initialize an EBC object.

        Args:
            matrix: a instance of SparseMatrix that represents the original distribution
            n_clusters: a list of number of clusters along each dimension
            max_iterations: the maximum number of iterations before we stop, default to be 10
            jitter_max: a small amount of value to add to probabilities to break ties when doing cluster assignments, default to be 1e-10
            objective_tolerance: the threshold of difference between two successive objective values for us to stop

        Return:
            None
        """
        if not isinstance(matrix, SparseMatrix):
            raise Exception("Matrix argument to EBC is not instance of SparseMatrix.")

        # check to ensure matrix is a probability distribution
        np.testing.assert_approx_equal(matrix.sum(), 1.0, significant=7,
                                       err_msg='Matrix elements does not sum to 1. Please normalize your matrix.')

        self.pXY = matrix  # the joint probability distribution e.g. p(X,Y)- the original sparse, multidimensional matrix
        self.K = n_clusters  # numbers of clusters along each dimension (len(K) = dim)
        self.dim = self.pXY.dim  # overall dimension of the matrix
        self.max_it = max_iterations

        self.cXY = None  # list of list: cluster assignments along each dimension e.g. [C(X), C(Y), ...]
        self.pX = None  # marginal probabilities

        self.qXhatYhat = None  # the approximate probability distribution after clustering e.g. q(X',Y'); need to be SparseMatrix
        self.qXhat = None  # the marginal clustering distribution in a list e.g. [q(X'), q(Y'), ...]
        self.qXxHat = None  # the distribution conditioned on the clustering in a list e.g. [q(X|X'), q(Y|Y'), ...]

        self.jitter_max = jitter_max  # amount to add to cluster assignment scores to break ties
        self.objective_tolerance = objective_tolerance  # the threshold for us to stop

    def run(self, assigned_clusters=None, verbose=True):
        """ To run the ebc algorithm.

        Args:
            assigned_clusters: an optional list of list representing the initial assignment of clusters along each dimension.

        Return:
            cXY: a list of list of cluster assignments along each dimension e.g. [C(X), C(Y), ...]
            objective: the final objective value
            max_it: the number of iterations that the algorithm has run
        """
        if verbose: print "Running EBC on a %d-d sparse matrix with size %s ..." % (self.dim, str(self.pXY.N))
        # Step 1: initialization steps
        self.pX = self.calculate_marginals(self.pXY)
        if assigned_clusters:
            if verbose: print "Using specified clusters, with cluster number on each axis: %s ..." % self.K
            self.cXY = assigned_clusters
        else:
            if verbose: print "Randomly initializing clusters, with cluster number on each axis: %s ..." % self.K
            self.cXY = self.initialize_cluster_centers(self.pXY, self.K)

        # Step 2: calculate cluster joint and marginal distributions
        self.qXhatYhat = self.calculate_joint_cluster_distribution(self.cXY, self.K, self.pXY)
        self.qXhat = self.calculate_marginals(self.qXhatYhat)  # the cluster marginals along each axis
        self.qXxHat = self.calculate_conditionals(self.cXY, self.pXY.N, self.pX, self.qXhat)

        # Step 3: iterate through dimensions, recalculating distributions
        last_objective = objective = INFINITE
        for t in xrange(self.max_it):
            if verbose: sys.stdout.write("--> Running iteration %d " % (t + 1)); sys.stdout.flush()
            for axis in xrange(self.dim):
                self.cXY[axis] = self.compute_clusters(self.pXY, self.qXhatYhat, self.qXhat, self.qXxHat, self.cXY,
                                                       axis)
                self.ensure_correct_number_clusters(self.cXY[axis], self.K[axis])  # check to ensure correct K
                self.qXhatYhat = self.calculate_joint_cluster_distribution(self.cXY, self.K, self.pXY)
                self.qXhat = self.calculate_marginals(self.qXhatYhat)
                self.qXxHat = self.calculate_conditionals(self.cXY, self.pXY.N, self.pX, self.qXhat)
                if verbose: sys.stdout.write("."); sys.stdout.flush()
            objective = self.calculate_objective()
            if verbose: sys.stdout.write(" objective value = %f\n" % (objective))
            if abs(objective - last_objective) < self.objective_tolerance:
                if verbose: print "EBC finished in %d iterations, with final objective value %.4f" % (t + 1, objective)
                return self.cXY, objective, t + 1
            last_objective = objective
        if verbose: print "EBC finished in %d iterations, with final objective value %.4f" % (self.max_it, objective)
        return self.cXY, objective, self.max_it  # hit max iterations - just return current assignments

    def compute_clusters(self, pXY, qXhatYhat, qXhat, qXxhat, cXY, axis):
        """ Compute the best cluster assignment along a single axis, given all the distributions and clusters on other axes.

        Args:
            pXY: the original probability distribution matrix
            qXhatYhat: the joint distribution over the clusters
            qXhat: the marginal distributions of qXhatYhat
            qXxhat: the distribution conditioned on the clustering in a list
            cXY: current cluster assignments along each dimension
            axis: the axis (dimension) over which clusters are being computed

        Return:
            Best cluster assignment along a single axis as a list
        """
        if not isinstance(pXY, SparseMatrix) or not isinstance(qXhatYhat, SparseMatrix):
            raise Exception("Arguments to compute_clusters not an instance of SparseMatrix.")
        # To assign clusters, we calculate argmin_xhat D(p(Y,Z|x) || q(Y,Z|xhat)),
        # where D(P|Q) = \sum_i P_i log (P_i / Q_i)
        dPQ = np.zeros(shape=(pXY.N[axis], qXhatYhat.N[axis]))
        # iterate though all non-zero elements; here we are making use of the sparsity to reduce computation
        for coords, p_i in pXY.nonzero_elements.iteritems():
            coord_this_axis = coords[axis]
            px = self.pX[axis][coord_this_axis]
            p_i = 1 if px == 0 else p_i / px  # calculate p(y|x) = p(x,y)/p(x), but we should be careful if px == 0
            current_cluster_assignments = [cXY[i][coords[i]] for i in
                                           xrange(self.dim)]  # cluster assignments on each axis
            for xhat in xrange(self.K[axis]):
                current_cluster_assignments[axis] = xhat  # temporarily assign dth dimension to this xhat
                current_qXhatYhat = qXhatYhat.get(tuple(current_cluster_assignments))
                current_qXhat = qXhat[axis][xhat]
                q_i = 1.0
                if current_qXhatYhat == 0 and current_qXhat == 0:
                    q_i = 0  # Here we define 0/0=0
                else:
                    q_i *= current_qXhatYhat / current_qXhat
                    for i in xrange(self.dim):
                        if i == axis: continue
                        q_i *= qXxhat[i][coords[i]]
                if q_i == 0:  # this can definitely happen if cluster joint distribution has zero element
                    dPQ[coord_this_axis, xhat] = INFINITE
                else:
                    dPQ[coord_this_axis, xhat] += p_i * math.log(p_i / q_i)

        # add random jitter to break ties
        dPQ += self.jitter_max * random_sample(dPQ.shape)
        return list(dPQ.argmin(1))  # return the closest cluster assignment under KL-divergence

    def calculate_marginals(self, pXY):
        """ Calculate the marginal probabilities given a joint distribution.

        Args:
            pXY: sparse matrix over which marginals are calculated, e.g. P(X,Y)

        Return:
            marginals: the calculated marginal probabilities, e.g. P(X)
        """
        if not isinstance(pXY, SparseMatrix):
            raise Exception("Illegal argument to marginal calculation: " + str(pXY))
        marginals = [[0] * Ni for Ni in pXY.N]
        for d in pXY.nonzero_elements:
            for i in xrange(len(d)):
                marginals[i][d[i]] += pXY.nonzero_elements[d]
        return marginals

    def calculate_joint_cluster_distribution(self, cXY, K, pXY):
        """ Calculate the joint cluster distribution q(X',Y') using the current prob distribution and
        cluster assignments. (Here we use X' to denote X_hat)

        Args:
            cXY: current cluster assignments for each axis
            K: numbers of clusters along each axis
            pXY: original probability distribution matrix

        Return:
            qXhatYhat: the joint cluster distribution
        """
        if not isinstance(pXY, SparseMatrix):
            raise Exception("Matrix argument to calculate_joint_cluster_distribution not an instance of SparseMatrix.")
        qXhatYhat = SparseMatrix(K)  # joint distribution over clusters
        for coords in pXY.nonzero_elements:
            # find the coordinates of the cluster for this element
            cluster_coords = []
            for i in xrange(len(coords)):
                cluster_coords.append(cXY[i][coords[i]])
            qXhatYhat.add_value(tuple(cluster_coords), pXY.nonzero_elements[coords])
        return qXhatYhat

    def calculate_conditionals(self, cXY, N, pX, qXhat):
        """ Calculate the conditional marginal distributions given the clustering distribution, i.e. q(X|X').

        Args:
            cXY: current cluster assignments
            N: lengths of each dimension in the original data matrix
            pX: marginal distributions over original data matrix
            qXhat: marginal distributions over cluster joint distribution

        Return:
            conditional_distributions: a list of distribution for each axis, with each element being a list of prob for i-th row/column in this axis.
        """
        conditional_distributions = [[0] * Ni for Ni in N]
        for i in xrange(len(cXY)):
            cluster_assignments_this_dimension = cXY[i]
            for j in xrange(len(cluster_assignments_this_dimension)):
                cluster = cluster_assignments_this_dimension[j]
                if pX[i][j] == 0 and qXhat[i][cluster] == 0:
                    conditional_distributions[i][j] = 0
                else:
                    conditional_distributions[i][j] = pX[i][j] / qXhat[i][cluster]
        return conditional_distributions

    def initialize_cluster_centers(self, pXY, K):
        """ Initializes the cluster assignments along each axis, by first selecting k centers, 
        and then map each row to its closet center under cosine similarity.

        Args:
            pXY: original data matrix
            K: numbers of clusters desired in each dimension

        Return:
            new_C: a list of list of cluster id that the current index in the current axis is assigned to.
        """
        if not isinstance(pXY, SparseMatrix):
            raise Exception("Matrix argument to initialize_cluster_centers is not an instance of SparseMatrix.")
        new_C = [[-1] * Ni for Ni in pXY.N]

        for axis in xrange(len(K)):  # loop over each dimension
            # choose cluster centers
            axis_length = pXY.N[axis]
            center_indices = random.sample(xrange(axis_length), K[axis])
            cluster_ids = {}
            for i in xrange(K[axis]):  # assign identifiers to clusters
                center_index = center_indices[i]
                cluster_ids[center_index] = i
            centers = defaultdict(lambda: defaultdict(float))  # all nonzero indices for each center
            for coords in pXY.nonzero_elements:
                coord_this_axis = coords[axis]
                if coord_this_axis in cluster_ids:  # is a center
                    reduced_coords = tuple(
                        [coords[i] for i in xrange(len(coords)) if i != axis])  # coords without the current axis
                    centers[cluster_ids[coord_this_axis]][reduced_coords] = pXY.nonzero_elements[
                        coords]  # (cluster_id, other coords) -> value

            # assign rows to clusters
            scores = np.zeros(shape=(pXY.N[axis], K[axis]))  # scores: axis_size x cluster_number
            denoms_P = np.zeros(shape=(pXY.N[axis]))
            denoms_Q = np.zeros(shape=(K[axis]))
            for coords in pXY.nonzero_elements:
                coord_this_axis = coords[axis]
                if coord_this_axis in center_indices:
                    continue  # don't reassign cluster centers, please
                reduced_coords = tuple([coords[i] for i in xrange(len(coords)) if i != axis])
                for cluster_index in cluster_ids:
                    xhat = cluster_ids[cluster_index]  # need cluster ID, not the axis index
                    if reduced_coords in centers[xhat]:  # overlapping point
                        P_i = pXY.nonzero_elements[coords]
                        Q_i = centers[xhat][reduced_coords]
                        scores[coords[axis]][xhat] += P_i * Q_i  # now doing based on cosine similarity
                        denoms_P[coords[axis]] += P_i * P_i  # magnitude of this slice of original matrix
                        denoms_Q[xhat] += Q_i * Q_i  # magnitude of cluster centers

            # normalize scores
            scores = divide(scores, outer(sqrt(denoms_P), sqrt(denoms_Q)))
            scores[scores == 0] = -1.0

            # add random jitter to scores to handle tie-breaking
            scores += self.jitter_max * random_sample(scores.shape)
            new_cXYi = list(scores.argmax(1))  # this needs to be argmax because cosine similarity

            # make sure to assign the cluster centers to themselves
            for center_index in cluster_ids:
                new_cXYi[center_index] = cluster_ids[center_index]

            # ensure numbers of clusters are correct
            self.ensure_correct_number_clusters(new_cXYi, K[axis])
            new_C[axis] = new_cXYi
        return new_C

    def ensure_correct_number_clusters(self, cXYi, expected_K):
        """ To ensure a cluster assignment actually has the expected total number of clusters. 

        Args:
            cXYi: the input cluster assignment
            expected_K: expected number of clusters on this axis

        Return:
            None. The assignment will be changed in place in cXYi.
        """
        clusters_represented = set()
        for c in cXYi:
            clusters_represented.add(c)
        if len(clusters_represented) == expected_K:
            return
        for c in xrange(expected_K):
            if c not in clusters_represented:
                index_to_change = random.randint(0, len(cXYi) - 1)
                cXYi[index_to_change] = c
        self.ensure_correct_number_clusters(cXYi, expected_K)

    def calculate_objective(self):
        """ Calculate the value of the objective function given the current cluster assignments. 

        Return:
            objective: the objective function value
        """
        objective = 0.0
        for d in self.pXY.nonzero_elements:
            pXY_i = self.pXY.nonzero_elements[d]
            qXY_i = self.get_element_approximate_distribution(d)
            if qXY_i == 0:
                print(pXY_i)
            objective += pXY_i * math.log(pXY_i / qXY_i)
        return objective

    def get_element_approximate_distribution(self, coords):
        """ Get the distribution approximated by q(X,Y). """
        clusters = [self.cXY[i][coords[i]] for i in xrange(len(coords))]
        element = self.qXhatYhat.get(tuple(clusters))
        for i in xrange(len(coords)):
            element *= self.qXxHat[i][coords[i]]
        return element


def main():
    """ An example run of EBC. """
    with open("resources/matrix-ebc-paper-sparse.tsv", "r") as f:
        data = []
        for line in f:
            sl = line.split("\t")
            if len(sl) < 5:  # headers
                continue
            data.append([sl[0], sl[2], float(sl[4])])

    matrix = SparseMatrix([14052, 7272])
    matrix.read_data(data)
    matrix.normalize()
    ebc = EBC(matrix, [30, 125], 10, 1e-10, 0.01)
    cXY, objective, it = ebc.run()


if __name__ == "__main__":
    main()
