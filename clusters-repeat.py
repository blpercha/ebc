from collections import defaultdict
import sys

from ebc import EBC
from matrix import SparseMatrix


def main():
    data_file = sys.argv[1]
    ebc_cols = [int(e) for e in sys.argv[2].split(",")]
    K = [int(e) for e in sys.argv[3].split(",")]
    N_runs = int(sys.argv[4])
    output_file = sys.argv[5]
    jitter_max = float(sys.argv[6])
    max_iterations_ebc = int(sys.argv[7])
    entity_cols = [int(e) for e in sys.argv[8].split(",")]
    object_toler = float(sys.argv[9])

    # get original data
    raw_data = [line.split("\t") for line in open(data_file, "r")]
    data = [[d[i] for i in ebc_cols] for d in raw_data]
    data_dimensions = len(data[0]) - 1

    # get axis length for each dimension
    N = []
    for dim in range(data_dimensions):
        N.append(len(set([d[dim] for d in data])))
    print(N)

    # set up matrix
    M = SparseMatrix(N)
    M.read_data(data)
    M.normalize()

    # set up entity map to ids
    entity_map = defaultdict(tuple)
    for d in raw_data:
        entity = tuple([d[i] for i in entity_cols])
        entity_ids = tuple([M.feature_ids[ebc_cols.index(i)][d[i]] for i in entity_cols])
        entity_map[entity_ids] = entity

    # figure out which ebc columns the entity columns correspond to
    entity_column_indices = []
    for c in ebc_cols:
        if c in entity_cols:
            entity_column_indices.append(ebc_cols.index(c))

    # run EBC and get entity cluster assignments
    ebc_M = EBC(M, K, max_iterations_ebc, jitter_max, object_toler)
    clusters = defaultdict(list)
    for t in range(N_runs):
        print "run ", t
        cXY_M, objective_M, it_M = ebc_M.run()
        for e1 in entity_map.keys():
            c1_i = tuple([cXY_M[i][e1[i]] for i in entity_column_indices])
            clusters[e1].append(c1_i)

    # print assignments
    writer = open(output_file, "w")
    for k in clusters:
        e1_name = entity_map[k]
        writer.write(",".join([str(e) for e in k]) + "\t" +
                     ",".join([e for e in e1_name]) + "\t" + "\t".join([",".join([str(f) for f in e])
                                                                        for e in clusters[k]]) + "\n")
        writer.flush()
    writer.close()


if __name__ == "__main__":
    main()
