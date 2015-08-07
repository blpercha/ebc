from collections import defaultdict
import sys

from numpy import zeros, arange
from numpy.linalg import matrix_rank
from scipy import linalg

data_file = sys.argv[1]
cols = [int(e) for e in sys.argv[2].split(",")]

coords = defaultdict(float)
xs = set()
ys = set()

f = open(data_file, "r")
for line in f:
    sl = line.split("\t")
    x = sl[cols[0]]
    y = sl[cols[1]]
    xs.add(x)
    ys.add(y)
    val = float(sl[cols[2]])
    coords[(x, y)] = val

f.close()

print len(xs), len(ys)

matrix = zeros(shape=[len(xs), len(ys)], dtype=float)

xs = list(xs)
ys = list(ys)
xs.sort()
ys.sort()
for i in range(len(xs)):
    for j in range(len(ys)):
        matrix[i, j] = coords[(xs[i], ys[j])]

U, s, Vh = linalg.svd(matrix)
for t in arange(0, 10, 0.05):
    row_rank = matrix_rank(matrix, tol=t)
    col_rank = matrix_rank(matrix.transpose(), tol=t)
    print t, row_rank, col_rank
