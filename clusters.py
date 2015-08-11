from collections import defaultdict
from copy import copy
from operator import itemgetter
from random import randint, random
import sys
from math import floor

from ebc import EBC
from matrix import SparseMatrix

data_file = sys.argv[1]
cols = [int(e) for e in sys.argv[2].split(",")]
N_iterations = int(sys.argv[3])
step_size = int(sys.argv[4])
n_avg = int(sys.argv[5])
output_file = open(sys.argv[6], "w")
dim_divisor = float(sys.argv[7])

raw_data = [line.split("\t") for line in open(data_file, "r")]
data = [[d[i] for i in cols] for d in raw_data]
data_dimensions = len(data[0]) - 1

# get axis length for each dimension
N = []
for dim in range(data_dimensions):
    N.append(len(set([d[dim] for d in data])))
print(N)


# optimize numbers of clusters for this matrix

def roundDown(x, to_nearest):
    return int(floor(x / float(to_nearest))) * to_nearest


K_old = [roundDown(randint(step_size, int(N[i] / dim_divisor)), step_size)
         for i in range(data_dimensions)]
D_0 = 1e10  # old objective difference


def pickNewK(k_old, max_dim, step, n, states_tried, axis_divisor):
    k = copy(k_old)
    while tuple(k) in states_tried:
        current_dim = randint(0, max_dim - 1)
        max_clusters = int(n[current_dim] / axis_divisor)
        k[current_dim] = k[current_dim] + step if random() < 0.5 else k[current_dim] - step
        if k[current_dim] <= 0:
            k[current_dim] = step
        elif k[current_dim] >= max_clusters:
            k[current_dim] = max_clusters
    return k


state_best = defaultdict(float)
max_states = 1
for N_i in N:
    max_states *= int(1.0 * N_i / (dim_divisor * step_size))
print "Max states: ", max_states

for j in range(N_iterations):
    if len(state_best) == max_states:
        break

    K = pickNewK(K_old, data_dimensions, step_size, N, state_best, dim_divisor)

    output_file.write(str(j) + "\t" + ",".join([str(e) for e in K]) + "\t" +
                      ",".join([str(e) for e in K_old]) + "\t")
    output_file.flush()
    print j, K, K_old,

    # because the data are so noisy, we average over many trials to get our estimate
    D_1 = 0.0
    for l in range(n_avg):
        M = SparseMatrix(N)
        M.read_data(data)
        Mr = M.shuffle_old()  # could also be M.shuffle()

        M.normalize()

        ebc_M = EBC(M, K, 100)
        cXY_M, objective_M = ebc_M.run()

        Mr.normalize()

        ebc_Mr = EBC(Mr, K, 100)
        cXY_Mr, objective_Mr = ebc_Mr.run()

        D_1_l = objective_M - objective_Mr
        D_1 += D_1_l
    D_1 /= n_avg

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
