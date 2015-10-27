from collections import defaultdict
import random

import numpy as np
from scipy.sparse import dok_matrix

class EBC2D:
    def __init__(self, matrix, n_clusters, max_iterations=10, jitter_max=1e-10, objective_tolerance=0.01):
        if not isinstance(matrix, dok_matrix):
            raise Exception("Matrix argument to EBC2D needs to be dok_matrix.")

        np.testing.assert_approx_equal(matrix.sum(), 1.0, significant=7, err_msg= \
            'Matrix elements does not sum to 1. Please normalize your matrix.')

        matrix = matrix.tocsr() # convert the matrix from dok_matrix to csc_matrix for speeding up
        self.pXY = matrix # p(X,Y)
        self.N = self.pXY.shape
        self.cX = np.zeros(self.N[0], dtype=np.int) # cluster assignment along X axis, C(X)
        self.cY = np.zeros(self.N[1], dtype=np.int) # C(Y)
        self.K = n_clusters # number of clusters along x and y axis
        self.max_iters = max_iterations

        self.pX = np.empty(self.N[0]) # marginal probabilities, p(X)
        self.pY = np.empty(self.N[1]) # marginal probabilities, p(Y)

        self.qXhatYhat = dok_matrix(tuple(self.K), dtype=np.float32) # the approx probability distribution after clustering
        self.qXhat = np.zeros(self.K[0]) # q(X')
        self.qYhat = np.zeros(self.K[1]) # q(Y')
        self.qX_xhat = np.zeros(self.N[0]) # q(X|X')
        self.qY_yhat = np.zeros(self.N[1]) # q(Y|Y')
        
        self.jitter_max = jitter_max
        self.objective_tolerance = objective_tolerance

    def run(self, assigned_clusters=None):
        # Step 1: initialization steps
        self.pX, self.pY = self.calculate_marginals(self.pXY)
        if assigned_clusters and len(assigned_clusters) == 2:
            self.cX = np.asarray(assigned_clusters[0])
            self.cY = np.asarray(assigned_clusters[1])
        else:
            self.cX, self.cY = self.initialize_cluster_centers(self.pXY, self.K)

        # print self.cX
        # print self.cY
        # Step 2: calculate cluster joint and marginal distributions

        # Step 3: iterate, recalculating distributions

        objective = 0

        return [self.cX, self.cY], objective, self.max_iters

    def calculate_marginals(self, pXY):
        """ Calculate the marginal probabilities given a joint distribution.

        Args:
            pXY: sparse [multidimensional] matrix over which marginals are calculated.
        """
        pX = pXY.sum(1) # sum along the y dimension, note that the dimension index should be reverse
        pY = pXY.sum(0) # sum along the x dimension
        return pX, pY

    def calcualate_conditionals(self, cX, cY, N, pX, pY, qXhat, qYhat):
        """ Calculate the conditional marginal distributions given the clustering distribution, i.e. q(X|X').

        Args:
            cX, cY: current cluster assignments
            N: lengths of each dimension in the original data matrix
            pX, pY: marginal distributions over original data matrix
            qXhat, qYhat: marginal distributions over cluster joint distribution

        Return:
            qX_xhat, qY_yhat: conditional marginal distributions for each axis.
        """
        qX_xhat = pX / qXhat[cX] # qX_xhat is a N[0]-size vector here
        qY_yhat = pY / qYhat[cY] # note that it could be problematic if cX is not a int vector
        # TODO: check for divide-by-zero errors
        return qX_xhat, qY_yhat

    def initialize_cluster_centers(self, pXY, K):
        """ Some magic code that initializes the cluster along each axis.

        Args:
            pXY: original data matrix
            K: numbers of clusters desired in each dimension

        Return:
            cX, cY: a list of cluster id that the current index in the current axis is assigned to.
        """
        # For x axis
        centers = pXY[random.sample(range(pXY.shape[0]), K[0]), :].toarray() # randomly select clustering centers
        cX = self.assign_clusters(pXY, centers, axis=0)
        # For y axis
        centers = pXY[:, random.sample(range(pXY.shape[1]), K[1])].toarray() # randomly select clustering centers
        cY = self.assign_clusters(pXY, centers, axis=1)
        # TODO: do I need to check to ensure the correct number of clusters?
        return cX, cY

    def assign_clusters(self, pXY, centers, axis):
        # TODO: fully vectorize the code
        scores = np.zeros(shape=(pXY.shape[axis], centers.shape[axis]))
        for i in range(pXY.shape[axis]):
            row = pXY[i, :].toarray() if axis == 0 else pXY[:, i].toarray() # here I need to force row to be np array to perform the deduction
            score_i = row * centers
            # TODO: any way to automatically broadcast a sparse csr matrix?
            score_i = score_i.sum(1) if axis == 0 else score_i.sum(0) # calculate u.v
            centers_length = np.power(centers, 2).sum(1) if axis == 0 else np.power(centers, 2).sum(0) # get |v|
            score_i /= centers_length # get u.v / |v|
            scores[i,:] = score_i.flatten()
        scores += self.jitter_max * np.random.mtrand.random_sample(scores.shape)
        C = np.argmax(scores, 1)
        return C

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
    m = dok_matrix((nrow, ncol), dtype=np.float32)
    for d in data:
        r = feature_ids[0][d[0]]
        c = feature_ids[1][d[1]]
        value = float(d[2])
        if value != 0.0:
            m[r, c] = value
    # normalize the matrix
    m = m / m.sum()
    return m


