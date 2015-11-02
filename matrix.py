from collections import defaultdict
from operator import itemgetter
from random import shuffle


class SparseMatrix:
    """ An implementation of sparse matrix that is used by the ITCC and EBC algorithm. """

    def __init__(self, N):
        """ Initialize the sparse matrix.

        Args:
            N: the size of the matrix on each axis in a list-like data structure
        """
        self.dim = len(N)  # dimensionality of matrix
        self.nonzero_elements = {}
        self.N = N
        # feature_ids should be a map from feature name to the corresponding index. 
        # For example, in a 2D matrix, each feature corresponds to a specific row or column.
        self.feature_ids = defaultdict(lambda: defaultdict(int))

    def read_data(self, data):
        """ Read the data from a list and populate the matrix. If 'data' is not a list, simply return. 
        
        Args:
            data: each element of the data list should be a list, and should have the following form:
            [feature1, feature2, ..., feature dim, value]
        """
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
        """ Get an element of the sparse matrix.

        Args:
            coordinates: indices of the element as a tuple

        Return:
            the element value
        """
        if coordinates in self.nonzero_elements:
            return self.nonzero_elements[coordinates]
        return 0.0

    def set(self, coordinates, value):
        """ Set the value for an element in the sparse matrix.

        Args:
            coordinates: indices of the element as a tuple
            value: the element value
        """
        self.nonzero_elements[coordinates] = value

    def add_value(self, coordinates, added_value):
        """ Add a specific value to an element in the sparse matrix.

        Args:
            coordinates: indices of the element as a tuple
            added_value: the value to add
        """
        if coordinates in self.nonzero_elements:
            self.nonzero_elements[coordinates] += added_value
        else:
            self.nonzero_elements[coordinates] = added_value

    def sum(self):
        """ Get the sum of all sparse matrix elements. 

        Return:
            the sum value
        """
        sum_values = 0.0
        for v in self.nonzero_elements.values():
            sum_values += v
        return sum_values

    def normalize(self):
        """ Normalize the sparse matrix such that the elements in the matrix sum up to 1. """
        sum_values = self.sum()
        for d in self.nonzero_elements:
            self.nonzero_elements[d] /= sum_values

    def shuffle(self):
        """ Randomly shuffle the nonzero elements in the original matrix,  and return a new matrix with the elements shuffled.

        Return:
            a new sparse matrix with all the elements shuffled
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

    def __str__(self):
        value_list = sorted(self.nonzero_elements.items(), key=itemgetter(0), reverse=False)
        return "\n".join(["\t".join([str(e) for e in v[0]]) + "\t" + str(v[1]) for v in value_list])
