import unittest

import numpy as np

from matrix import SparseMatrix
from ebc import EBC
import ebc2d
from ebc2d import EBC2D


class TestSanityCheck(unittest.TestCase):
    """ Do a sanity check for the EBC code, using the data from the original ITCC paper. """

    def setUp(self):
        with open("resources/matrix-itcc-paper-orig.tsv", "r") as f:
            data = [l.split('\t') for l in f]

        self.matrix = SparseMatrix([6, 6])
        self.matrix.read_data(data)
        self.matrix.normalize()

    def cartesian(self, arrays, out=None):
        arrays = [np.asarray(x) for x in arrays]
        dtype = arrays[0].dtype

        n = np.prod([x.size for x in arrays])
        if out is None:
            out = np.zeros([n, len(arrays)], dtype=dtype)

        m = n / arrays[0].size
        out[:, 0] = np.repeat(arrays[0], m)
        if arrays[1:]:
            self.cartesian(arrays[1:], out=out[0:m, 1:])
            for j in xrange(1, arrays[0].size):
                out[j * m:(j + 1) * m, 1:] = out[0:m, 1:]
        return out

    def testEbcOnSparseMatrix(self):
        ebc = EBC(self.matrix, [3, 2], 10, 1e-10, 0.01)
        cXY, objective, it = ebc.run(verbose=False)
        print "--> ebc"
        print "objective: ", objective
        print "iterations: ", it

        ebc = EBC(self.matrix, [3, 2], 10, 1e-10, 0.01)
        ebc.run(assigned_clusters=[[2, 0, 1, 1, 2, 2], [0, 0, 1, 0, 1, 1]], verbose=False)
        indices = [range(N_d) for N_d in ebc.pXY.N]
        index_list = self.cartesian(indices)
        approx_distribution = {}
        for location in index_list:
            q = 1.0
            c_location = []
            for i in range(len(location)):
                c_i = ebc.cXY[i][location[i]]
                c_location.append(c_i)
                q *= ebc.qXxHat[i][location[i]]
            q *= ebc.qXhatYhat.get(tuple(c_location))
            approx_distribution[tuple(location)] = q

        self.assertAlmostEquals(approx_distribution[(0, 0)], 0.054)
        self.assertAlmostEquals(approx_distribution[(0, 1)], 0.054)
        self.assertAlmostEquals(approx_distribution[(0, 2)], 0.042)
        self.assertAlmostEquals(approx_distribution[(0, 3)], 0.0)
        self.assertAlmostEquals(approx_distribution[(0, 4)], 0.0)
        self.assertAlmostEquals(approx_distribution[(0, 5)], 0.0)
        self.assertAlmostEquals(approx_distribution[(1, 0)], 0.054)
        self.assertAlmostEquals(approx_distribution[(1, 1)], 0.054)
        self.assertAlmostEquals(approx_distribution[(1, 2)], 0.042)
        self.assertAlmostEquals(approx_distribution[(1, 3)], 0.0)
        self.assertAlmostEquals(approx_distribution[(1, 4)], 0.0)
        self.assertAlmostEquals(approx_distribution[(1, 5)], 0.0)
        self.assertAlmostEquals(approx_distribution[(2, 0)], 0.0)
        self.assertAlmostEquals(approx_distribution[(2, 1)], 0.0)
        self.assertAlmostEquals(approx_distribution[(2, 2)], 0.0)
        self.assertAlmostEquals(approx_distribution[(2, 3)], 0.042)
        self.assertAlmostEquals(approx_distribution[(2, 4)], 0.054)
        self.assertAlmostEquals(approx_distribution[(2, 5)], 0.054)
        self.assertAlmostEquals(approx_distribution[(3, 0)], 0.0)
        self.assertAlmostEquals(approx_distribution[(3, 1)], 0.0)
        self.assertAlmostEquals(approx_distribution[(3, 2)], 0.0)
        self.assertAlmostEquals(approx_distribution[(3, 3)], 0.042)
        self.assertAlmostEquals(approx_distribution[(3, 4)], 0.054)
        self.assertAlmostEquals(approx_distribution[(3, 5)], 0.054)
        self.assertAlmostEquals(approx_distribution[(4, 0)], 0.036)
        self.assertAlmostEquals(approx_distribution[(4, 1)], 0.036)
        self.assertAlmostEquals(approx_distribution[(4, 2)], 0.028)
        self.assertAlmostEquals(approx_distribution[(4, 3)], 0.028)
        self.assertAlmostEquals(approx_distribution[(4, 4)], 0.036)
        self.assertAlmostEquals(approx_distribution[(4, 5)], 0.036)
        self.assertAlmostEquals(approx_distribution[(5, 0)], 0.036)
        self.assertAlmostEquals(approx_distribution[(5, 1)], 0.036)
        self.assertAlmostEquals(approx_distribution[(5, 2)], 0.028)
        self.assertAlmostEquals(approx_distribution[(5, 3)], 0.028)
        self.assertAlmostEquals(approx_distribution[(5, 4)], 0.036)
        self.assertAlmostEquals(approx_distribution[(5, 5)], 0.036)

    def testEbc2dOnSparseMatrix(self):
        with open("resources/matrix-itcc-paper-orig.tsv", "r") as f:
            data = [l.split('\t') for l in f]
        m = ebc2d.get_matrix_from_data(data)
        # run without assigned clusters
        ebc = EBC2D(m, [3, 2], 10, 1e-10, 0.01)
        cXY, objective, it = ebc.run(verbose=False)
        print "--> ebc2d"
        print "objective: ", objective
        print "iterations: ", it

        # run with assigned clusters
        ebc = EBC2D(m, [3, 2], 10, 1e-10, 0.01)
        cXY, objective, it = ebc.run(assigned_clusters=[[2, 0, 1, 1, 2, 2], [0, 0, 1, 0, 1, 1]], verbose=False)
        indices = [range(N_d) for N_d in ebc.pXY.shape]
        index_list = self.cartesian(indices)
        approx_distribution = {}
        qX_xhat = [ebc.qX_xhat, ebc.qY_yhat]
        for location in index_list:
            q = 1.0
            c_location = []
            for i in range(len(location)):
                c_i = cXY[i][location[i]]
                c_location.append(c_i)
                q *= qX_xhat[i][location[i]]
            q *= ebc.qXhatYhat[c_location[0], c_location[1]]
            approx_distribution[tuple(location)] = q

        self.assertAlmostEquals(approx_distribution[(0, 0)], 0.054)
        self.assertAlmostEquals(approx_distribution[(0, 1)], 0.054)
        self.assertAlmostEquals(approx_distribution[(0, 2)], 0.042)
        self.assertAlmostEquals(approx_distribution[(0, 3)], 0.0)
        self.assertAlmostEquals(approx_distribution[(0, 4)], 0.0)
        self.assertAlmostEquals(approx_distribution[(0, 5)], 0.0)
        self.assertAlmostEquals(approx_distribution[(1, 0)], 0.054)
        self.assertAlmostEquals(approx_distribution[(1, 1)], 0.054)
        self.assertAlmostEquals(approx_distribution[(1, 2)], 0.042)
        self.assertAlmostEquals(approx_distribution[(1, 3)], 0.0)
        self.assertAlmostEquals(approx_distribution[(1, 4)], 0.0)
        self.assertAlmostEquals(approx_distribution[(1, 5)], 0.0)
        self.assertAlmostEquals(approx_distribution[(2, 0)], 0.0)
        self.assertAlmostEquals(approx_distribution[(2, 1)], 0.0)
        self.assertAlmostEquals(approx_distribution[(2, 2)], 0.0)
        self.assertAlmostEquals(approx_distribution[(2, 3)], 0.042)
        self.assertAlmostEquals(approx_distribution[(2, 4)], 0.054)
        self.assertAlmostEquals(approx_distribution[(2, 5)], 0.054)
        self.assertAlmostEquals(approx_distribution[(3, 0)], 0.0)
        self.assertAlmostEquals(approx_distribution[(3, 1)], 0.0)
        self.assertAlmostEquals(approx_distribution[(3, 2)], 0.0)
        self.assertAlmostEquals(approx_distribution[(3, 3)], 0.042)
        self.assertAlmostEquals(approx_distribution[(3, 4)], 0.054)
        self.assertAlmostEquals(approx_distribution[(3, 5)], 0.054)
        self.assertAlmostEquals(approx_distribution[(4, 0)], 0.036)
        self.assertAlmostEquals(approx_distribution[(4, 1)], 0.036)
        self.assertAlmostEquals(approx_distribution[(4, 2)], 0.028)
        self.assertAlmostEquals(approx_distribution[(4, 3)], 0.028)
        self.assertAlmostEquals(approx_distribution[(4, 4)], 0.036)
        self.assertAlmostEquals(approx_distribution[(4, 5)], 0.036)
        self.assertAlmostEquals(approx_distribution[(5, 0)], 0.036)
        self.assertAlmostEquals(approx_distribution[(5, 1)], 0.036)
        self.assertAlmostEquals(approx_distribution[(5, 2)], 0.028)
        self.assertAlmostEquals(approx_distribution[(5, 3)], 0.028)
        self.assertAlmostEquals(approx_distribution[(5, 4)], 0.036)
        self.assertAlmostEquals(approx_distribution[(5, 5)], 0.036)
