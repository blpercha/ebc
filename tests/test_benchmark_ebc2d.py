import unittest

import ebc2d
from ebc2d import EBC2D


class TestBenchmarkEBC2D(unittest.TestCase):
    """ Benchmark the EBC2D code as a unittest, using the sparse matrix data. """

    def setUp(self):
        with open("resources/matrix-ebc-paper-sparse.tsv", "r") as f:
            data = []
            for line in f:
                sl = line.split("\t")
                if len(sl) < 5:  # headers
                    continue
                data.append([sl[0], sl[2], float(sl[4])])

        self.matrix = ebc2d.get_matrix_from_data(data)

    def testEbc2dOnSparseMatrix(self):
        ebc = EBC2D(self.matrix, [30, 125], 10, 1e-10, 0.01)
        cXY, objective, it = ebc.run()
        print "objective: ", objective
        print "iterations: ", it
        # self.assertEquals(len(ebc.pXY.nonzero[0]), 29456)
        self.assertEquals(len(set(cXY[0])), 30)
        self.assertEquals(len(set(cXY[1])), 125)
