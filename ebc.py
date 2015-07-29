from math import log

from numpy import zeros, asarray, prod, repeat
from scipy.lib.six import xrange

from matrix import SparseMatrix


class EBC:
    # matrix: input SparseMatrix
    # K: numbers of clusters along each dimension (length is number of dimensions)
    def __init__(self, matrix, K, max_iterations):
        if not isinstance(matrix, SparseMatrix):
            return
        self.pXY = matrix  # the original sparse, multidimensional matrix
        self.cXY = [[0] * Ni for Ni in self.pXY.N]  # cluster assignments each dimension
        self.K = K  # numbers of clusters along each dimension (len(K) = D)
        self.D = self.pXY.D  # overall dimension of the matrix
        self.max_iterations = max_iterations

        self.pXY_marginals = [[0] * Ni for Ni in self.pXY.N]

        self.qXhatYhat = SparseMatrix(self.D, K)
        self.qXhatYhat_marginals = [[0] * Ki for Ki in self.K]
        self.qXxHat = [[0] * Ni for Ni in self.pXY.N]

    def run_ebc(self, assigned_C=None):
        # Step 1: initialization steps
        self.pXY_marginals = self.calculate_marginals(self.pXY.M, self.pXY.N)
        self.cXY = self.initialize_cluster_centers(self.pXY.M, self.pXY.N, assigned_C)

        # Step 2: calculate cluster joint and marginal distributions
        self.qXhatYhat = self.calculate_joint_cluster_distribution(self.cXY, self.K, self.D, self.pXY.M)
        self.qXhatYhat_marginals = self.calculate_marginals(self.qXhatYhat.M, self.K)
        self.qXxHat = self.calculate_conditionals(self.cXY, self.pXY.N, self.pXY_marginals, self.qXhatYhat_marginals)

        # Step 3: iterate through dimensions, recalculating distributions
        last_cXY = self.cXY.copy()
        for t in range(self.max_iterations):
            for d in range(self.D):
                self.cXY[d] = self.compute_clusters(self.pXY, self.qXhatYhat, self.qXhatYhat_marginals,
                                                    self.qXxHat, self.cXY, d)
                self.qXhatYhat = self.calculate_joint_cluster_distribution(self.cXY, self.K, self.D, self.pXY.M)
                self.qXhatYhat_marginals = self.calculate_marginals(self.qXhatYhat.M, self.K)
                self.qXxHat = self.calculate_conditionals(self.cXY, self.pXY.N, self.pXY_marginals,
                                                          self.qXhatYhat_marginals)
            if self.cXY == last_cXY:
                break
            else:
                print(self.cXY)
                last_cXY = self.cXY.copy()

    def get_approximate_distribution(self):
        indices = [range(N_d) for N_d in self.pXY.N]
        index_list = self.cartesian(indices)
        approx_distribution = {}
        for location in index_list:
            q = 1.0
            c_location = []
            for i in range(len(location)):
                c_i = self.cXY[i][location[i]]
                c_location.append(c_i)
                q *= self.qXxHat[i][location[i]]
            q *= self.qXhatYhat.get(tuple(c_location))
            approx_distribution[tuple(location)] = q
        return approx_distribution

    def cartesian(self, arrays, out=None):
        arrays = [asarray(x) for x in arrays]
        dtype = arrays[0].dtype

        n = prod([x.size for x in arrays])
        if out is None:
            out = zeros([n, len(arrays)], dtype=dtype)

        m = n / arrays[0].size
        out[:, 0] = repeat(arrays[0], m)
        if arrays[1:]:
            self.cartesian(arrays[1:], out=out[0:m, 1:])
            for j in xrange(1, arrays[0].size):
                out[j * m:(j + 1) * m, 1:] = out[0:m, 1:]
        return out

    """
    :param sparse_P: the original data matrix
    :param sparse_Q: the joint distribution over the clusters
    :param marginals_Q: the marginal distributions of sparse_Q
    :param cond_P: for each dimension, the conditional distributions over the clusters (in a single vector)
    :param C: current cluster assignments
    :param d: the dimension over which clusters are being computed
    """

    def compute_clusters(self, sparse_P, sparse_Q, marginals_Q, cond_P, C, d):
        # want argmax_xhat D(p(Y,Z|x) || q(Y,Z|xhat))
        # D(P|Q) = \sum_i P_i log (P_i / Q_i)
        dPQ = zeros(shape=(sparse_P.N[d], sparse_Q.N[d]))
        for coords in sparse_P.M:
            x = coords[d]
            clust = [C[i][coords[i]] for i in range(len(coords))]
            P_i = sparse_P.M[coords]
            for xhat in range(sparse_Q.N[d]):
                clust[d] = xhat  # temporarily assign dth dimension to this xhat
                Q_i = 1.0
                for i in range(len(clust)):
                    if i == d:
                        continue
                    Q_i *= cond_P[i][coords[i]]
                Q_i *= sparse_Q.get(tuple(clust)) / marginals_Q[d][xhat]
                if Q_i == 0:  # this can definitely happen if cluster joint distribution has zero element
                    dPQ[x, xhat] = 1e10
                else:
                    dPQ[x, xhat] += P_i * log(P_i / Q_i)
        return list(dPQ.argmin(1))

    def calculate_marginals(self, sparse_M, N):
        marginals = [[0] * Ni for Ni in N]
        for d in sparse_M:
            for i in range(len(d)):
                marginals[i][d[i]] += sparse_M[d]
        return marginals

    def calculate_joint_cluster_distribution(self, C, K, dim, sparse_M):
        cluster_joint_distribution = SparseMatrix(dim, K)
        for d in sparse_M:
            c = []
            for i in range(len(d)):
                c.append(C[i][d[i]])
            cluster_joint_distribution.add(tuple(c), sparse_M[d])
        return cluster_joint_distribution

    def calculate_conditionals(self, C, N, p_marg, c_marg):
        conditional_distributions = [[0] * Ni for Ni in N]
        for i in range(len(C)):
            cluster_assignments_this_dimension = C[i]
            for j in range(len(cluster_assignments_this_dimension)):
                cluster = cluster_assignments_this_dimension[j]
                conditional_distributions[i][j] = p_marg[i][j] / c_marg[i][cluster]
        return conditional_distributions

    def calculate_crosses(self):
        cross_distributions = [SparseMatrix(2, [self.pXY.N[i], self.qXhatYhat.N[i]]) for i in range(len(self.K))]
        for c in cross_distributions:
            print(c.to_string(), end="\n")

    def initialize_cluster_centers(self, sparse_matrix, N, assigned_C=None):
        if assigned_C:
            return assigned_C
        new_C = [[0] * Ni for Ni in self.pXY.N]
        return new_C  # TODO: ensure number of clusters is correct
