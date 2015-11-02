from collections import defaultdict
import random
import sys

import numpy as np


class EBC2D:
    def __init__(self, matrix, n_clusters, max_iterations=10, jitter_max=1e-10, objective_tolerance=0.01):
        """ To initialize an EBC object.

        Args:
            matrix: a instance of numpy ndarray that represents the original distribution
            n_clusters: a list of number of clusters along each dimension
            max_iterations: the maximum number of iterations before we stop, default to be 10
            jitter_max: a small amount of value to add to probabilities to break ties when doing cluster assignments, default to be 1e-10
            objective_tolerance: the threshold of difference between two successive objective values for us to stop

        Return:
            None
        """
        if not isinstance(matrix, np.ndarray):
            raise Exception("Matrix argument to EBC2D needs to be a numpy ndarray.")

        np.testing.assert_approx_equal(matrix.sum(), 1.0, significant=3, err_msg= \
            'Matrix elements does not sum to 1. Please normalize your matrix.')

        self.pXY = matrix  # p(X,Y)
        self.N = self.pXY.shape
        self.cX = np.zeros(self.N[0], dtype=np.int)  # cluster assignment along X axis, C(X)
        self.cY = np.zeros(self.N[1], dtype=np.int)  # C(Y)
        self.K = n_clusters  # number of clusters along x and y axis
        self.max_iters = max_iterations

        self.pX = np.empty(self.N[0])  # marginal probabilities, p(X)
        self.pY = np.empty(self.N[1])  # marginal probabilities, p(Y)

        self.qXhatYhat = np.zeros(tuple(self.K),
                                  dtype=np.float32)  # the approx probability distribution after clustering
        self.qXhat = np.zeros(self.K[0])  # q(X')
        self.qYhat = np.zeros(self.K[1])  # q(Y')
        self.qX_xhat = np.zeros(self.N[0])  # q(X|X')
        self.qY_yhat = np.zeros(self.N[1])  # q(Y|Y')

        self.jitter_max = jitter_max
        self.objective_tolerance = objective_tolerance

    def run(self, assigned_clusters=None, verbose=True):
        """ To run the ebc algorithm.

        Args:
            assigned_clusters: an optional list of list representing the initial assignment of clusters along each dimension.

        Return:
            cXY: a list of lists of cluster assignments along each dimension e.g. [C(X), C(Y), ...]
            objective: the final objective value
            max_it: the number of iterations that the algorithm has run
        """
        if verbose: print "Running EBC2D on a 2-d matrix with size %s ..." % str(self.pXY.shape)
        # Step 1: initialization steps
        self.pX, self.pY = self.calculate_marginals(self.pXY)
        if assigned_clusters and len(assigned_clusters) == 2:
            if verbose: print "Using specified clusters, with cluster number on each axis: %s ..." % self.K
            self.cX = np.asarray(assigned_clusters[0])
            self.cY = np.asarray(assigned_clusters[1])
        else:
            if verbose: print "Randomly initializing clusters, with cluster number on each axis: %s ..." % self.K
            self.cX, self.cY = self.initialize_cluster_centers(self.pXY, self.K)

        # Step 2: calculate cluster joint and marginal distributions
        self.qXhatYhat = self.calculate_joint_cluster_distribution(self.cX, self.cY, self.K, self.pXY)
        self.qXhat, self.qYhat = self.calculate_marginals(self.qXhatYhat)
        self.qX_xhat, self.qY_yhat = self.calculate_conditionals(self.cX, self.cY, self.pX, self.pY, self.qXhat,
                                                                 self.qYhat)

        # Step 3: iterate, recalculating distributions
        last_objective = objective = 1e10
        for it in range(self.max_iters):
            if verbose: sys.stdout.write("--> Running iteration %d " % (it + 1)); sys.stdout.flush()
            # compute row/column clusters
            for axis in range(2):
                if axis == 0:
                    self.cX = self.compute_row_clusters(self.pXY, self.qXhatYhat, self.qXhat, self.qY_yhat, self.cY)
                else:
                    self.cY = self.compute_col_clusters(self.pXY, self.qXhatYhat, self.qYhat, self.qX_xhat, self.cX)
                if verbose: sys.stdout.write("."); sys.stdout.flush()
            self.qXhatYhat = self.calculate_joint_cluster_distribution(self.cX, self.cY, self.K, self.pXY)
            self.qXhat, self.qYhat = self.calculate_marginals(self.qXhatYhat)
            self.qX_xhat, self.qY_yhat = self.calculate_conditionals(self.cX, self.cY, self.pX, self.pY, self.qXhat,
                                                                     self.qYhat)

            objective = self.calculate_objective()
            if verbose: sys.stdout.write(" objective value = %f\n" % (objective))
            if abs(objective - last_objective) < self.objective_tolerance:
                if verbose: print "EBC2D finished in %d iterations, with final objective value %.4f" % (
                it + 1, objective)
                return [self.cX.tolist(), self.cY.tolist()], objective, it + 1
            last_objective = objective
        if verbose: print "EBC2D finished in %d iterations, with final objective value %.4f" % (
        self.max_iters, objective)
        return [self.cX.tolist(), self.cY.tolist()], objective, self.max_iters

    def compute_row_clusters(self, pXY, qXhatYhat, qXhat, qY_yhat, cY):
        """ Compute the best row cluster assignment, given all the distributions and clusters on y axis.

        Args:
            pXY: the original probability distribution matrix
            qXhatYhat: the joint distribution over the clusters
            qXhat: the marginal distributions of qXhatYhat
            qY_yhat: the distribution conditioned on the y axis clustering in a list
            cY: current cluster assignments along y axis

        Return:
            Best row cluster assignment as a list
        """
        nrow, nc_row = pXY.shape[0], qXhat.shape[0]
        dPQ = np.empty((nrow, nc_row))
        # Step 1: generate q(y|x'): |x'| x |y|
        # - first expand q(x'y') to |x'| x |y| using clustering information cY
        expanded_qXhatYhat = np.nan_to_num((qXhatYhat.T / qXhat).T[:, cY])
        qY_xhat = expanded_qXhatYhat * qY_yhat
        # Step 2: loop through all clusters
        for i in range(nc_row):
            for j in range(nrow):
                pXY_row = pXY[j, :]
                with np.errstate(divide='ignore', invalid='ignore'):
                    log_matrix = np.log(pXY_row / qY_xhat[i, :])
                    log_matrix[log_matrix == -np.inf] = 0
                    log_matrix = np.nan_to_num(log_matrix)
                    dPQ[j, i] = pXY_row.dot(log_matrix.T)

        dPQ += self.jitter_max * np.random.mtrand.random_sample(dPQ.shape)
        C = dPQ.argmin(1)
        self.ensure_correct_number_clusters(C, nc_row)
        return C

    def compute_col_clusters(self, pXY, qXhatYhat, qYhat, qX_xhat, cX):
        """ Compute the best column cluster assignment, given all the distributions and clusters on x axis. """
        ncol, nc_col = pXY.shape[1], qYhat.shape[0]
        dPQ = np.empty((nc_col, ncol))
        expanded_qXhatYhat = np.nan_to_num((qXhatYhat / qYhat)[cX, :])
        qX_yhat = expanded_qXhatYhat.T * qX_xhat
        for i in range(nc_col):
            for j in range(ncol):
                pXY_col = pXY[:, j].T
                with np.errstate(divide='ignore', invalid='ignore'):
                    log_matrix = np.log(pXY_col / qX_yhat[i, :])
                    log_matrix[log_matrix == -np.inf] = 0
                    log_matrix = np.nan_to_num(log_matrix)
                    dPQ[i, j] = pXY_col.dot(log_matrix.T)
        dPQ += self.jitter_max * np.random.mtrand.random_sample(dPQ.shape)
        C = dPQ.argmin(0)
        self.ensure_correct_number_clusters(C, nc_col)
        return C

    def calculate_marginals(self, pXY):
        """ Calculate the marginal probabilities given a joint distribution.

        Args:
            pXY: distribution over which marginals are calculated.
        """
        pX = pXY.sum(1)  # sum along the y dimension, note that the dimension index should be reverse
        pY = pXY.sum(0)  # sum along the x dimension
        return np.squeeze(np.asarray(pX)), np.squeeze(np.asarray(pY))  # return a numpy array

    def calculate_conditionals(self, cX, cY, pX, pY, qXhat, qYhat):
        """ Calculate the conditional marginal distributions given the clustering distribution, i.e. q(X|X').

        Args:
            cX, cY: current cluster assignments
            pX, pY: marginal distributions over original data matrix
            qXhat, qYhat: marginal distributions over cluster joint distribution

        Return:
            qX_xhat, qY_yhat: conditional marginal distributions for each axis.
        """
        with np.errstate(divide='ignore', invalid='ignore'):
            qX_xhat = pX / qXhat[cX]  # qX_xhat is a N[0]-size vector here
            qY_yhat = pY / qYhat[cY]  # note that it could be problematic if cX is not a int vector
            qX_xhat[qX_xhat == np.inf] = 0
            qY_yhat[qY_yhat == np.inf] = 0
        return qX_xhat, qY_yhat  # want a Nx1 array-like matrix

    def calculate_joint_cluster_distribution(self, cX, cY, K, pXY):
        """ Calculate the joint cluster distribution q(X',Y') = p(X',Y') using the current prob distribution and
        cluster assignments. (Here we use X' to denote X_hat)

        Args:
            cX, cY: current cluster assignments for each axis
            K: numbers of clusters along each axis
            pXY: original probability distribution matrix

        Return:
            qXhatYhat: the joint cluster distribution
        """
        nc_row, nc_col = K  # num of clusters along row and col
        qXhatYhat = np.zeros(K)
        # itm_matrix = csr_matrix((nc_row, pXY.shape[1])) # nc_row * col sparse intermidiate matrix
        itm_matrix = np.empty((nc_row, pXY.shape[1]))
        for i in range(nc_row):
            itm_matrix[i, :] = pXY[np.where(cX == i)[0], :].sum(0)
        for i in range(nc_col):
            qXhatYhat[:, i] = itm_matrix[:, np.where(cY == i)[0]].sum(1).flatten()
        return qXhatYhat

    def initialize_cluster_centers(self, pXY, K):
        """ Initializes the cluster assignments along each axis, by first selecting k centers, 
        and then map each row to its closet center under cosine similarity.

        Args:
            pXY: original data matrix
            K: numbers of clusters desired in each dimension

        Return:
            cX, cY: a list of cluster id that the current index in the current axis is assigned to.
        """
        # For x axis
        centers = pXY[random.sample(range(pXY.shape[0]), K[0]), :]  # randomly select clustering centers
        cX = self.assign_clusters(pXY, centers, axis=0)
        self.ensure_correct_number_clusters(cX, K[0])
        # For y axis
        centers = pXY[:, random.sample(range(pXY.shape[1]), K[1])]  # randomly select clustering centers
        cY = self.assign_clusters(pXY, centers, axis=1)
        self.ensure_correct_number_clusters(cY, K[1])
        return cX, cY  # return a numpy array

    def assign_clusters(self, pXY, centers, axis):
        """ Assign each row/col to clusters given cluster centers on this axis. """
        cluster_num = centers.shape[axis]
        scores = np.zeros(shape=(pXY.shape[axis], cluster_num))
        for i in range(cluster_num):
            if axis == 0:
                center = centers[i, :]
                score_i = pXY * center
            else:
                center = centers[:, i]
                score_i = pXY.T * center
            score_i = score_i.sum(1)  # calculate u.v
            center_length = np.linalg.norm(center)  # get |v|
            score_i = score_i / center_length  # get u.v/|v|
            scores[:, i] = score_i.flatten()
        scores += self.jitter_max * np.random.mtrand.random_sample(scores.shape)
        C = np.argmax(scores, 1)
        return C

    def calculate_objective(self):
        """ Calculate the KL-divergence between p(X,Y) and q(X,Y). 
        Here q(x,y) can be written as p(x',y')*p(x|x')*p(y|y'). """
        # Here I cannot vectorize the computation
        objective = .0
        x_indices, y_indices = np.nonzero(self.pXY)
        # compute values for all useful elements in qXY
        for i in range(len(x_indices)):
            x_idx, y_idx = x_indices[i], y_indices[i]
            v = self.pXY[x_idx, y_idx]
            c_x, c_y = self.cX[x_idx], self.cY[y_idx]
            v_qXY = self.qX_xhat[x_idx] * self.qY_yhat[y_idx] * self.qXhatYhat[c_x, c_y]
            objective += v * np.log(v / v_qXY)
        return objective

    def ensure_correct_number_clusters(self, C, expected_K):
        """ To ensure a cluster assignment actually has the expected total number of clusters. 

        Args:
            cXYi: the input cluster assignment
            expected_K: expected number of clusters on this axis

        Return:
            None. The assignment will be changed in place.
        """
        clusters_unique = np.unique(C)
        num_clusters = clusters_unique.shape[0]
        if num_clusters == expected_K:
            return
        for c in range(expected_K):
            if num_clusters < c + 1 or clusters_unique[c] != c:  # no element assigned to c
                idx = random.randint(0, C.shape[0] - 1)
                C[idx] = c
        self.ensure_correct_number_clusters(C, expected_K)


def get_matrix_from_data(data):
    """ Read the data from a list and construct a scipy sparse dok_matrix. If 'data' is not a list, simply return. 
    
    Each element of the data list should be a list, and should have the following form:
        [feature1, feature2, ..., feature dim, value]
    """
    feature_ids = defaultdict(lambda: defaultdict(int))
    for d in data:
        location = []
        for i in range(len(d) - 1):
            f_i = d[i]
            if f_i not in feature_ids[i]:
                feature_ids[i][f_i] = len(feature_ids[i])  # new index is size of dict
            location.append(feature_ids[i][f_i])
    nrow = len(feature_ids[0])
    ncol = len(feature_ids[1])
    m = np.zeros((nrow, ncol), dtype=np.float32)
    for d in data:
        r = feature_ids[0][d[0]]
        c = feature_ids[1][d[1]]
        value = float(d[2])
        if value != 0.0:
            m[r, c] = value
    # normalize the matrix
    m = m / m.sum()
    return m
