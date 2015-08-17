from collections import defaultdict
from math import log
import random
from copy import copy

from numpy import zeros
from numpy.ma import divide, outer, sqrt
from numpy.random.mtrand import random_sample
from numpy.testing import assert_approx_equal

from matrix import SparseMatrix


class EBC:
    def __init__(self, matrix, n_clusters, max_iterations, jitter_max):
        if not isinstance(matrix, SparseMatrix):
            raise Exception("Matrix argument to EBC is not SparseMatrix.")

        # check to ensure matrix is a probability distribution
        sum_values = 0.0
        for nonzero_value in matrix.nonzero_elements.values():
            sum_values += nonzero_value
        assert_approx_equal(sum_values, 1.0, significant=7)

        self.pXY = matrix  # the original sparse, multidimensional matrix
        self.cXY = [[0] * Ni for Ni in self.pXY.N]  # cluster assignments each dimension
        self.K = n_clusters  # numbers of clusters along each dimension (len(K) = D)
        self.dim = self.pXY.dim  # overall dimension of the matrix
        self.max_it = max_iterations

        self.pX = None

        self.qXhatYhat = SparseMatrix(self.K)
        self.qXhat = [[0] * Ki for Ki in self.K]
        self.qXxHat = [[0] * Ni for Ni in self.pXY.N]

        # amount to add to cluster assignment scores to break ties
        self.jitter_max = jitter_max

    def run(self, assigned_C=None):
        # Step 1: initialization steps
        self.pX = self.calculate_marginals(self.pXY)
        self.cXY = self.initialize_cluster_centers(self.pXY, self.K, assigned_C)

        # Step 2: calculate cluster joint and marginal distributions
        self.qXhatYhat = self.calculate_joint_cluster_distribution(self.cXY, self.K, self.pXY)
        self.qXhat = self.calculate_marginals(self.qXhatYhat)
        self.qXxHat = self.calculate_conditionals(self.cXY, self.pXY.N, self.pX, self.qXhat)

        # Step 3: iterate through dimensions, recalculating distributions
        last_cXY = copy(self.cXY)
        for t in range(self.max_it):
            K_order = [d for d in range(self.dim)]
            random.shuffle(K_order)

            for dim in K_order:
                self.cXY[dim] = self.compute_clusters(self.pXY, self.qXhatYhat, self.qXhat, self.qXxHat, self.cXY, dim)
                self.ensure_correct_number_clusters(self.cXY[dim], self.K[dim])  # check to ensure correct K
                self.qXhatYhat = self.calculate_joint_cluster_distribution(self.cXY, self.K, self.pXY)
                self.qXhat = self.calculate_marginals(self.qXhatYhat)
                self.qXxHat = self.calculate_conditionals(self.cXY, self.pXY.N, self.pX, self.qXhat)
            if self.cXY == last_cXY:
                return self.cXY, self.calculate_objective()
            last_cXY = copy(self.cXY)
        return self.cXY, self.calculate_objective()  # hit max iterations - just return current assignments

    """
    :param pXY: the original data matrix
    :param qXhatYhat: the joint distribution over the clusters
    :param qXhat: the marginal distributions of qXhatYhat
    :param qXxhat: for each dimension, the conditional distributions over the clusters (in a single vector)
    :param cXY: current cluster assignments
    :param axis: the axis (dimension) over which clusters are being computed
    """

    def compute_clusters(self, pXY, qXhatYhat, qXhat, qXxhat, cXY, axis):
        if not isinstance(pXY, SparseMatrix) or not isinstance(qXhatYhat, SparseMatrix):
            raise Exception("Arguments to compute_clusters not sparse.")
        # want argmin_xhat D(p(Y,Z|x) || q(Y,Z|xhat))
        # D(P|Q) = \sum_i P_i log (P_i / Q_i)
        dPQ = zeros(shape=(pXY.N[axis], qXhatYhat.N[axis]))
        for coords in pXY.nonzero_elements:
            coordinate_this_axis = coords[axis]
            current_cluster_assignments = [cXY[i][coords[i]] for i in range(len(coords))]
            P_i = pXY.nonzero_elements[coords]
            for xhat in range(qXhatYhat.N[axis]):
                current_cluster_assignments[axis] = xhat  # temporarily assign dth dimension to this xhat
                Q_i = 1.0
                for i in range(len(current_cluster_assignments)):
                    if i == axis:
                        continue
                    Q_i *= qXxhat[i][coords[i]]
                if qXhatYhat.get(tuple(current_cluster_assignments)) == 0 and qXhat[axis][xhat] == 0:
                    Q_i = 0
                else:
                    Q_i *= qXhatYhat.get(tuple(current_cluster_assignments)) / qXhat[axis][xhat]
                if Q_i == 0:  # this can definitely happen if cluster joint distribution has zero element
                    dPQ[coordinate_this_axis, xhat] = 1e10
                else:
                    dPQ[coordinate_this_axis, xhat] += P_i * log(P_i / Q_i)

        # add random jitter to break ties
        dPQ += self.jitter_max * random_sample(dPQ.shape)
        return list(dPQ.argmin(1))

    """
    :param pXY: sparse [multidimensional] matrix over which marginals are calculated
    """

    def calculate_marginals(self, pXY):
        if not isinstance(pXY, SparseMatrix):
            raise Exception("Illegal argument to marginal calculation: " + str(pXY))
        marginals = [[0] * Ni for Ni in pXY.N]
        for d in pXY.nonzero_elements:
            for i in range(len(d)):
                marginals[i][d[i]] += pXY.nonzero_elements[d]
        return marginals

    """
    :param cXY: current cluster assignments
    :param K: numbers of clusters along each axis
    :param pXY: original data matrix
    """

    def calculate_joint_cluster_distribution(self, cXY, K, pXY):
        if not isinstance(pXY, SparseMatrix):
            raise Exception("Matrix argument to calculate_joint_cluster_distribution not sparse.")
        qXhatYhat = SparseMatrix(K)  # joint distribution over clusters
        for coords in pXY.nonzero_elements:
            cluster_coords = []
            for i in range(len(coords)):
                cluster_coords.append(cXY[i][coords[i]])
            qXhatYhat.add_value(tuple(cluster_coords), pXY.nonzero_elements[coords])
        return qXhatYhat

    """
    :param cXY: current cluster assignments
    :param N: lengths of each dimension in the original data matrix
    :param pX: marginal distributions over original data matrix
    :param qXhat: marginal distributions over cluster joint distribution
    """

    def calculate_conditionals(self, cXY, N, pX, qXhat):
        conditional_distributions = [[0] * Ni for Ni in N]
        for i in range(len(cXY)):
            cluster_assignments_this_dimension = cXY[i]
            for j in range(len(cluster_assignments_this_dimension)):
                cluster = cluster_assignments_this_dimension[j]
                if pX[i][j] == 0 and qXhat[i][cluster] == 0:
                    conditional_distributions[i][j] = 0
                else:
                    conditional_distributions[i][j] = pX[i][j] / qXhat[i][cluster]
        return conditional_distributions

    """
    :param pXY: original data matrix
    :param K: numbers of clusters desired in each dimension
    :param assigned_C: (optional) fixed initial assignments of clusters
    """

    def initialize_cluster_centers(self, pXY, K, assigned_C=None):
        if assigned_C:
            return assigned_C
        if not isinstance(pXY, SparseMatrix):
            raise Exception("Matrix argument to initialize_cluster_centers is not sparse.")
        new_C = [[-1] * Ni for Ni in pXY.N]

        for axis in range(len(K)):
            # choose cluster centers
            axis_length = pXY.N[axis]
            center_indices = random.sample(range(axis_length), K[axis])
            cluster_ids = {}
            for i in range(len(center_indices)):  # assign identifiers to clusters
                center_index = center_indices[i]
                cluster_ids[center_index] = i
            centers = defaultdict(lambda: defaultdict(float))  # all nonzero indices for each center
            for coords in pXY.nonzero_elements:
                index_on_axis = coords[axis]
                if index_on_axis in cluster_ids:
                    reduced_coords = tuple([coords[i] for i in range(len(coords)) if i != axis])
                    centers[cluster_ids[index_on_axis]][reduced_coords] = pXY.nonzero_elements[coords]

            # assign rows to clusters
            scores = zeros(shape=(pXY.N[axis], K[axis]))
            denoms_P = zeros(shape=(pXY.N[axis]))
            denoms_Q = zeros(shape=(K[axis]))
            for coords in pXY.nonzero_elements:
                index_on_axis = coords[axis]
                if index_on_axis in center_indices:
                    continue  # don't reassign cluster centers, please
                reduced_coords = tuple([coords[i] for i in range(len(coords)) if i != axis])
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
            new_cXYi = list(scores.argmax(1))  # !!!!!!! this needs to be argmax b/c cosine similarity!!!!

            # make sure to assign the cluster centers to themselves
            for center_index in cluster_ids:
                new_cXYi[center_index] = cluster_ids[center_index]

            # ensure numbers of clusters are correct
            self.ensure_correct_number_clusters(new_cXYi, K[axis])
            new_C[axis] = new_cXYi
        return new_C

    def ensure_correct_number_clusters(self, cXYi, expected_K):
        clusters_represented = set()
        for c in cXYi:
            clusters_represented.add(c)
        if len(clusters_represented) == expected_K:
            return
        for c in range(expected_K):
            if c not in clusters_represented:
                index_to_change = random.randint(0, len(cXYi) - 1)
                cXYi[index_to_change] = c
        self.ensure_correct_number_clusters(cXYi, expected_K)

    def calculate_objective(self):
        objective = 0.0
        for d in self.pXY.nonzero_elements:
            pXY_element = self.pXY.nonzero_elements[d]
            qXY_element = self.get_element_approximate_distribution(d)
            if qXY_element == 0:
                print(pXY_element)
            objective += pXY_element * log(pXY_element / qXY_element)
        return objective

    def get_element_approximate_distribution(self, coords):
        clusters = [self.cXY[i][coords[i]] for i in range(len(coords))]
        element = self.qXhatYhat.get(tuple(clusters))
        for i in range(len(coords)):
            element *= self.qXxHat[i][coords[i]]
        return element
