import unittest

from matrix import SparseMatrix


class TestMatrix(unittest.TestCase):
    def setUp(self):
        self.data = [l.split('\t') for l in open('sample-matrix-file.txt', 'r').readlines()]

    def testMatrixInit(self):
        matrix = SparseMatrix([2, 4, 9])
        matrix.read_data(self.data)
        self.assertEquals(matrix.nonzero_elements[(1, 3, 7)], 2.0)
        self.assertEquals(matrix.nonzero_elements[(0, 0, 0)], 2.0)
        self.assertEquals(matrix.nonzero_elements[(0, 0, 2)], 2.0)
        self.assertEquals(matrix.nonzero_elements[(1, 1, 5)], 7.0)
        self.assertEquals(matrix.nonzero_elements[(1, 1, 3)], 3.0)
        self.assertEquals(matrix.nonzero_elements[(1, 3, 6)], 2.0)
        self.assertEquals(matrix.nonzero_elements[(1, 3, 8)], 2.0)
        self.assertEquals(matrix.nonzero_elements[(0, 0, 1)], 2.0)
        self.assertEquals(matrix.nonzero_elements[(1, 1, 4)], 2.0)
        self.assertEquals(matrix.nonzero_elements[(1, 2, 5)], 2.0)
        self.assertEquals(len(matrix.nonzero_elements), 10)
        self.assertEquals(matrix.feature_ids[0], {'mice': 1, 'patient': 0})
        self.assertEquals(matrix.feature_ids[1], {'R92Q': 1, 'R91W': 2, 'Val30Met': 0, 'R90W': 3})
        self.assertEquals(matrix.feature_ids[2], {'START_ENTITY|nmod|END_ENTITY': 1,
                                                  'START_ENTITY|nummod|END_ENTITY': 5,
                                                  'FAP|compound|END_ENTITY': 2,
                                                  'expression|nmod|END_ENTITY': 8,
                                                  '+|compound|END_ENTITY': 7,
                                                  'mice|nummod|END_ENTITY': 3,
                                                  'homozygous|nsubj|START_ENTITY': 6,
                                                  'mutation|appos|END_ENTITY': 4,
                                                  'START_ENTITY|nmod|FAP': 0})
