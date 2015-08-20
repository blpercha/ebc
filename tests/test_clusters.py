import unittest
from ebc import EBC
from matrix import SparseMatrix


class TestMatrix(unittest.TestCase):
    def setUp(self):
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

        ebc = EBC(matrix, [3, 3, 3], 10, 1e-10)
        assigned_C = [[0, 0, 1, 1, 2, 2], [0, 0, 1, 1, 2, 2], [0, 0, 1, 1, 2, 2]]
        cXY, objective = ebc.run(assigned_C)
        self.assertEquals(cXY, assigned_C)
        self.assertAlmostEqual(objective, 0.0)
        cXY, objective = ebc.run()  # random initialization
        self.assertAlmostEqual(objective, 0.0)