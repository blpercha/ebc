from collections import defaultdict
from copy import copy
from operator import itemgetter
from random import randint, random
import sys

from ebc import EBC
from matrix import SparseMatrix

data_file = sys.argv[1]
cols = [int(e) for e in sys.argv[2].split(",")]
N_iterations = int(sys.argv[3])
step_size = int(sys.argv[4])
n_avg = int(sys.argv[5])
output_file = open(sys.argv[6], "w")

raw_data = [line.split("\t") for line in open(data_file, "r")]
data = [[d[i] for i in cols] for d in raw_data]
data_dimensions = len(data[0]) - 1

# get axis length for each dimension
N = []
for dim in range(data_dimensions):
    N.append(len(set([d[dim] for d in data])))
print(N)

# optimize numbers of clusters for this matrix
K_old = [randint(2, min(100, N[i])) for i in range(data_dimensions)]
D_0 = 1e10  # old objective difference

# keep map of best value for state
state_best = defaultdict(float)
for j in range(N_iterations):
    K = copy(K_old)
    dim_chosen = randint(0, data_dimensions - 1)
    if K[dim_chosen] >= N[dim_chosen] - step_size:  # at right boundary
        K[dim_chosen] -= step_size
    elif K[dim_chosen] <= step_size:  # at left boundary
        K[dim_chosen] += step_size
    else:
        K[dim_chosen] = min(K[dim_chosen] + step_size, N[dim_chosen]) \
            if random() < 0.5 else max(K[dim_chosen] - step_size, 2)

    output_file.write(str(j) + "\t" + ",".join([str(e) for e in K]) + "\t" +
                      ",".join([str(e) for e in K_old]) + "\t")
    output_file.flush()
    print j, K, K_old,

    # because the data are so noisy, we average over many trials to get our estimate
    D_1 = 0.0
    for l in range(n_avg):
        M = SparseMatrix(N)
        M.read_data(data)
        Mr = M.shuffle()  # could also be M.shuffle_disperse()

        M.normalize()

        ebc_M = EBC(M, K, 100)
        cXY_M, objective_M = ebc_M.run()

        Mr.normalize()

        ebc_Mr = EBC(Mr, K, 100)
        cXY_Mr, objective_Mr = ebc_Mr.run()

        D_1_l = objective_M - objective_Mr
        D_1 += D_1_l
    D_1 /= n_avg

    if D_1 < state_best[tuple(K)]:
        state_best[tuple(K)] = D_1

    output_file.write(str(D_1) + "\t" + str(D_0) + "\t")
    output_file.flush()
    print D_1, D_0,

    moved = False
    if D_1 < D_0:  # getting better
        K_old = K
        D_0 = D_1
        moved = True

    output_file.write(str(moved) + "\n")
    output_file.flush()
    print moved

for state_b in sorted(state_best.items(), key=itemgetter(1), reverse=False):
    print "\t".join([str(e) for e in state_b[0]]) + "\t" + str(state_b[1])
