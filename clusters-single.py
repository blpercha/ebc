import sys

from numpy import mean, std

from ebc import EBC
from matrix import SparseMatrix

data_file = sys.argv[1]
cols = [int(e) for e in sys.argv[2].split(",")]
K = [int(e) for e in sys.argv[3].split(",")]
N_trials = int(sys.argv[4])
output_file = sys.argv[5]
jitter_max = float(sys.argv[6])

# get original data
raw_data = [line.split("\t") for line in open(data_file, "r")]
data = [[d[i] for i in cols] for d in raw_data]
data_dimensions = len(data[0]) - 1

# get axis length for each dimension
N = []
for dim in range(data_dimensions):
    N.append(len(set([d[dim] for d in data])))
print(N)

D_1 = []
for j in range(N_trials):
    print "Trial ", j

    M = SparseMatrix(N)
    M.read_data(data)
    Mr = M.shuffle()  # could also be M.shuffle_old()

    M.normalize()

    ebc_M = EBC(M, K, 20, jitter_max)
    cXY_M, objective_M = ebc_M.run()

    Mr.normalize()

    ebc_Mr = EBC(Mr, K, 20, jitter_max)
    cXY_Mr, objective_Mr = ebc_Mr.run()

    D_1.append(objective_M - objective_Mr)

# write final result to combined file (other processes also write to this file)
output_stream = open(output_file, "a")
output_stream.write("\t".join([str(e) for e in K]) + "\t" + str(mean(D_1)) + "\t" + str(std(D_1)) + "\n")
output_stream.flush()
output_stream.close()
