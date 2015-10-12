from operator import itemgetter
import unittest

from numpy import asarray, prod, zeros, repeat

from scipy.lib.six import xrange

from ebc import EBC
from matrix import SparseMatrix


class TestEbc(unittest.TestCase):
    def setUp(self):
        self.data = [["0", "0", 0.05],
                     ["0", "1", 0.05],
                     ["0", "2", 0.05],
                     ["0", "3", 0.00],
                     ["0", "4", 0.00],
                     ["0", "5", 0.00],
                     ["1", "0", 0.05],
                     ["1", "1", 0.05],
                     ["1", "2", 0.05],
                     ["1", "3", 0.00],
                     ["1", "4", 0.00],
                     ["1", "5", 0.00],
                     ["2", "0", 0.00],
                     ["2", "1", 0.00],
                     ["2", "2", 0.00],
                     ["2", "3", 0.05],
                     ["2", "4", 0.05],
                     ["2", "5", 0.05],
                     ["3", "0", 0.00],
                     ["3", "1", 0.00],
                     ["3", "2", 0.00],
                     ["3", "3", 0.05],
                     ["3", "4", 0.05],
                     ["3", "5", 0.05],
                     ["4", "0", 0.04],
                     ["4", "1", 0.04],
                     ["4", "2", 0.00],
                     ["4", "3", 0.04],
                     ["4", "4", 0.04],
                     ["4", "5", 0.04],
                     ["5", "0", 0.04],
                     ["5", "1", 0.04],
                     ["5", "2", 0.04],
                     ["5", "3", 0.00],
                     ["5", "4", 0.04],
                     ["5", "5", 0.04]]
        self.matrix = SparseMatrix([6, 6])
        self.matrix.read_data(self.data)

    def testDataLoad(self):
        self.assertEquals(sorted(self.matrix.nonzero_elements.items(), key=itemgetter(0)),
                          [((0, 0), 0.05), ((0, 1), 0.05), ((0, 2), 0.05), ((1, 0), 0.05), ((1, 1), 0.05),
                           ((1, 2), 0.05), ((2, 3), 0.05), ((2, 4), 0.05), ((2, 5), 0.05), ((3, 3), 0.05),
                           ((3, 4), 0.05), ((3, 5), 0.05), ((4, 0), 0.04), ((4, 1), 0.04), ((4, 3), 0.04),
                           ((4, 4), 0.04), ((4, 5), 0.04), ((5, 0), 0.04), ((5, 1), 0.04), ((5, 2), 0.04),
                           ((5, 4), 0.04), ((5, 5), 0.04)])

    def testApproximateDistributionOriginalITCCPaper(self):
        ebc = EBC(self.matrix, [3, 2], 10, 1e-10, 0.01)
        ebc.run(assigned_C=[[2, 0, 1, 1, 2, 2], [0, 0, 1, 0, 1, 1]])
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

    def testOldMatrix(self):
        f = open("/Users/beth/Documents/phd/ebc/resources/matrix-ebc-paper-dense.tsv", "r")
        data = []
        for line in f:
            sl = line.split("\t")
            if len(sl) < 5:  # headers
                continue
            data.append([sl[0], sl[2], float(sl[4])])
        f.close()

        matrix = SparseMatrix([3514, 1232])
        matrix.read_data(data)
        matrix.normalize()
        ebc = EBC(matrix, [30, 125], 10, 1e-10, 0.01)
        cXY, objective, it = ebc.run()
        print "objective: ", objective
        print "iterations: ", it
        self.assertEquals(len(ebc.pXY.nonzero_elements), 10007)
        self.assertEquals(len(set(ebc.cXY[0])), 30)
        self.assertEquals(len(set(ebc.cXY[1])), 125)

    def testOldMatrix3d(self):
        f = open("/Users/beth/Documents/phd/ebc/resources/matrix-ebc-paper-dense-3d.tsv", "r")
        data = []
        for line in f:
            sl = line.split("\t")
            data.append([sl[0], sl[1], sl[2], float(sl[3])])
        f.close()

        matrix = SparseMatrix([756, 996, 1232])
        matrix.read_data(data)
        matrix.normalize()
        ebc = EBC(matrix, [30, 30, 10], 100, 1e-10, 0.01)
        cXY, objective, it = ebc.run()
        print "objective: ", objective
        print "iterations: ", it
        self.assertEquals(len(ebc.pXY.nonzero_elements), 10007)
        self.assertEquals(len(set(ebc.cXY[0])), 30)
        self.assertEquals(len(set(ebc.cXY[1])), 30)
        self.assertEquals(len(set(ebc.cXY[2])), 10)

    def test3DMatrix(self):
        data = [[0, 0, 0, 1.0],
                [0, 0, 1, 1.0],
                [0, 1, 0, 1.0],
                [0, 1, 1, 1.0],
                [1, 0, 0, 1.0],
                [1, 0, 1, 1.0],
                [1, 1, 0, 1.0],
                [1, 1, 1, 1.0],
                [2, 2, 2, 1.0],
                [2, 2, 3, 1.0],
                [2, 3, 2, 1.0],
                [3, 2, 2, 1.0],
                [2, 3, 3, 1.0],
                [3, 3, 2, 1.0],
                [3, 2, 3, 1.0],
                [3, 3, 3, 1.0],
                [4, 4, 4, 1.0],
                [4, 4, 5, 1.0],
                [4, 5, 4, 1.0],
                [4, 5, 5, 1.0],
                [5, 4, 4, 1.0],
                [5, 4, 5, 1.0],
                [5, 5, 4, 1.0],
                [5, 5, 5, 1.0]]
        matrix = SparseMatrix([6, 6, 6])
        matrix.read_data(data)
        matrix.normalize()
        ebc = EBC(matrix, [3, 3, 3], 10, 1e-10, 0.01)
        assigned_C = [[0, 0, 1, 1, 2, 2], [0, 0, 1, 1, 2, 2], [0, 0, 1, 1, 2, 2]]
        cXY, objective, it = ebc.run(assigned_C)
        self.assertEquals(cXY, assigned_C)
        self.assertAlmostEqual(objective, 0.0)
        self.assertEquals(it, 1)

        for i in range(100):
            cXY, objective, it = ebc.run()  # random initialization
            print cXY, objective, it
