from collections import defaultdict
from operator import itemgetter


class SparseMatrix:
    def __init__(self, dim, N):
        self.D = dim  # dimensionality of matrix
        self.M = {}
        self.N = N
        self.n_nonzero = 0
        self.feature_ids = defaultdict(lambda: defaultdict(int))

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
                self.M[tuple(location)] = value
        self.n_nonzero = len(self.M)

    def get(self, coordinates):
        if coordinates in self.M:
            return self.M[coordinates]
        return 0.0

    def set(self, coordinates, value):
        self.M[coordinates] = value

    def add(self, coordinates, added_value):
        if coordinates in self.M:
            self.M[coordinates] += added_value
        else:
            self.M[coordinates] = added_value

    def to_string(self):
        value_list = sorted(self.M.items(), key=itemgetter(0), reverse=False)
        return "\n".join(["\t".join([str(e) for e in v[0]]) + "\t" + str(v[1]) for v in value_list])
