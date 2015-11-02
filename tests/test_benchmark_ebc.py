import unittest

from matrix import SparseMatrix
from ebc import EBC


class TestBenchmarkEBC(unittest.TestCase):
    """ Benchmark the EBC code as a unittest, using the sparse matrix data. """

    def setUp(self):
        with open("resources/matrix-ebc-paper-sparse.tsv", "r") as f:
            data = []
            for line in f:
                sl = line.split("\t")
                if len(sl) < 5:  # headers
                    continue
                data.append([sl[0], sl[2], float(sl[4])])

        self.matrix = SparseMatrix([14052, 7272])
        self.matrix.read_data(data)
        self.matrix.normalize()

    def testEbcOnSparseMatrix(self):
        ebc = EBC(self.matrix, [30, 125], 10, 1e-10, 0.01)
        cXY, objective, it = ebc.run()
        print "objective: ", objective
        print "iterations: ", it
        self.assertEquals(len(ebc.pXY.nonzero_elements), 29456)
        self.assertEquals(len(set(ebc.cXY[0])), 30)
        self.assertEquals(len(set(ebc.cXY[1])), 125)
