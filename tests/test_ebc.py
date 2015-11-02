from operator import itemgetter
import unittest

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

    def testOldMatrix(self):
        with open("resources/matrix-ebc-paper-dense.tsv", "r") as f:
            data = []
            for line in f:
                sl = line.split("\t")
                if len(sl) < 5:  # headers
                    continue
                data.append([sl[0], sl[2], float(sl[4])])

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
        with open("resources/matrix-ebc-paper-dense-3d.tsv", "r") as f:
            data = []
            for line in f:
                sl = line.split("\t")
                data.append([sl[0], sl[1], sl[2], float(sl[3])])

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
