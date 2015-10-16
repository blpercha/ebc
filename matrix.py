from collections import defaultdict
from operator import itemgetter
from random import shuffle

class SparseMatrix:
    """ An implementation of sparse matrix that is used by the ITCC and EBC algorithm.
    """
    def __init__(self, N):
        self.dim = len(N)  # dimensionality of matrix
        self.nonzero_elements = {}
        self.N = N
        self.feature_ids = defaultdict(lambda: defaultdict(int))

    def read_data(self, data):
        """ Read the data from a list and populate the matrix. If 'data' is not a list, simply return. """
        if not isinstance(data, list):  # we expect a list of data points
            return
        for d in data:
            location = []
            for i in range(len(d) - 1):
                f_i = d[i]
                if f_i not in self.feature_ids[i]:
                    self.feature_ids[i][f_i] = len(self.feature_ids[i])  # new index is size of dict
                location.append(self.feature_ids[i][f_i])
            value = float(d[len(d) - 1])
            if value != 0.0:
                self.nonzero_elements[tuple(location)] = value

    def get(self, coordinates):
        if coordinates in self.nonzero_elements:
            return self.nonzero_elements[coordinates]
        return 0.0

    def set(self, coordinates, value):
        self.nonzero_elements[coordinates] = value

    def add_value(self, coordinates, added_value):
        if coordinates in self.nonzero_elements:
            self.nonzero_elements[coordinates] += added_value
        else:
            self.nonzero_elements[coordinates] = added_value

    def normalize(self):
        """ Normalize the sparse matrix such that the elements in the matrix sum up to 1. """
        sum_values = 0.0
        for d in self.nonzero_elements:
            sum_values += self.nonzero_elements[d]
        for d in self.nonzero_elements:
            self.nonzero_elements[d] /= sum_values

    def shuffle(self):
        """ Randomly shuffle the nonzero elements in the original matrix, 
        and return a new matrix with the elements shuffled.
        """
        self_shuffled = SparseMatrix(self.N)
        indices = []
        # Get all the indices of nonzero elements. indices is a list of 'dim' lists, each being a list of indices for a specific dimension
        for i in range(self.dim):
            indices.append([e[i] for e in self.nonzero_elements])
        for i in range(self.dim):
            shuffle(indices[i])
        values = [self.nonzero_elements[e] for e in self.nonzero_elements]
        shuffle(values)
        for j in range(len(self.nonzero_elements)):
            self_shuffled.add_value(tuple([indices[i][j] for i in range(self.dim)]), values[j])
        return self_shuffled

    def to_string(self):
        value_list = sorted(self.nonzero_elements.items(), key=itemgetter(0), reverse=False)
        return "\n".join(["\t".join([str(e) for e in v[0]]) + "\t" + str(v[1]) for v in value_list])
