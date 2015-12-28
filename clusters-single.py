import sys

from numpy import mean, std

from ebc import EBC
from matrix import SparseMatrix


def compareRandom(num_trials, tensor_dimensions, matrix_data, cluster_dimensions,
                  maxit_ebc, jitter_max_ebc, objective_tolerance):
    deltas = []
    objectives_M = []
    objectives_Mr = []
    iterations_M = []
    iterations_Mr = []
    noconverge_M = 0
    noconverge_Mr = 0
    for j in range(num_trials):
        print "Trial ", j

        M = SparseMatrix(tensor_dimensions)
        M.read_data(matrix_data)
        Mr = M.shuffle()  # could also be M.shuffle_old()

        M.normalize()

        ebc_M = EBC(M, cluster_dimensions, maxit_ebc, jitter_max_ebc, objective_tolerance)
        cXY_M, objective_M, it_M = ebc_M.run()
        if it_M == maxit_ebc:
            noconverge_M += 1
        else:
            iterations_M.append(it_M)

        Mr.normalize()

        ebc_Mr = EBC(Mr, cluster_dimensions, maxit_ebc, jitter_max_ebc, objective_tolerance)
        cXY_Mr, objective_Mr, it_Mr = ebc_Mr.run()
        if it_Mr == maxit_ebc:
            noconverge_Mr += 1
        else:
            iterations_Mr.append(it_Mr)

        objectives_M.append(objective_M)
        objectives_Mr.append(objective_Mr)
        deltas.append(objective_M - objective_Mr)
    return deltas, objectives_M, objectives_Mr, iterations_M, iterations_Mr, noconverge_M, noconverge_Mr


def main():
    data_file = sys.argv[1]
    cols = [int(e) for e in sys.argv[2].split(",")]
    K = [int(e) for e in sys.argv[3].split(",")]
    N_trials = int(sys.argv[4])
    output_file = sys.argv[5]
    jitter_max = float(sys.argv[6])
    max_iterations_ebc = int(sys.argv[7])
    object_tol = float(sys.argv[8])

    # get original data
    raw_data = [line.split("\t") for line in open(data_file, "r")]
    data = [[d[i] for i in cols] for d in raw_data]
    data_dimensions = len(data[0]) - 1

    # get axis length for each dimension
    N = []
    for dim in range(data_dimensions):
        N.append(len(set([d[dim] for d in data])))
    print(N)

    D_1, obj_orig, obj_rand, it_orig, it_rand, noconv_orig, noconv_rand = compareRandom(num_trials=N_trials,
                                                                                        tensor_dimensions=N,
                                                                                        matrix_data=data,
                                                                                        cluster_dimensions=K,
                                                                                        maxit_ebc=max_iterations_ebc,
                                                                                        jitter_max_ebc=jitter_max,
                                                                                        objective_tolerance=object_tol)

    # write final result to combined file (other processes also write to this file)
    output_stream = open(output_file, "a")
    output_stream.write("\t".join([str(e) for e in K]) + "\t" + str(mean(D_1)) + "\t" + str(std(D_1)) +
                        "\t" + str(mean(obj_orig)) + "\t" + str(mean(obj_rand)) +
                        "\t" + str(mean(it_orig)) + "\t" + str(mean(it_rand)) +
                        "\t" + str(noconv_orig) + "\t" + str(noconv_rand) + "\n")
    output_stream.flush()
    output_stream.close()


if __name__ == "__main__":
    main()
