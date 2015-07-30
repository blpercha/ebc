from collections import defaultdict
from operator import itemgetter


class SparseMatrix:
    def __init__(self, N):
        self.dim = len(N)  # dimensionality of matrix
        self.nonzero_elements = {}
        self.N = N
        self.feature_ids = defaultdict(lambda: defaultdict(int))
        self.sum_values = 0.0

    def read_data(self, data):
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
                self.sum_values += value

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
        for d in self.nonzero_elements:
            self.nonzero_elements[d] /= self.sum_values

    def to_string(self):
        value_list = sorted(self.nonzero_elements.items(), key=itemgetter(0), reverse=False)
        return "\n".join(["\t".join([str(e) for e in v[0]]) + "\t" + str(v[1]) for v in value_list])
