from random import randint
import sys

from ebc import EBC

from matrix import SparseMatrix

data_file = sys.argv[1]
cols = [int(e) for e in sys.argv[2].split(",")]
cluster_opt_iterations = int(sys.argv[3])

raw_data = [line.split("\t") for line in open(data_file, "r")]
data = [[d[i] for i in cols] for d in raw_data]
data_dimensions = len(data[0]) - 1

# get axis length for each dimension
N = []
for dim in range(data_dimensions):
    N.append(len(set([d[dim] for d in data])))

# optimize numbers of clusters for this matrix
K_old = [2 for i in range(data_dimensions)]
K = K_old
old_objective_difference = 1e10
for t in range(cluster_opt_iterations):
    print(t, K, old_objective_difference)

    M = SparseMatrix(N)
    M.read_data(data)
    Mr = M.shuffle()

    M.normalize()

    ebc_M = EBC(M, K, 100)
    cXY_M, objective_M = ebc_M.run()

    Mr.normalize()

    ebc_Mr = EBC(Mr, K, 100)
    cXY_Mr, objective_Mr = ebc_Mr.run()

    objective_difference = objective_M - objective_Mr
    if objective_difference < old_objective_difference:  # getting better
        K_old = K
        K = [max(2, K[i] + randint(-5, 5)) for i in range(data_dimensions)]
    else:
        K_old_copy = K_old
        K_old = K
        K = K_old_copy
    old_objective_difference = objective_difference
