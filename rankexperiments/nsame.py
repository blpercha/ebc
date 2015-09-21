from random import shuffle

from numpy import mean, array

from matrix import SparseMatrix

nrow = 14052
ncol = 7272

raw_data = [l.split('\t') for l in open('../resources/matrix-ebc-paper-sparse.tsv', 'r').readlines()]
data = [[e[0], e[2], e[4]] for e in raw_data]
matrix = SparseMatrix([nrow, ncol])
matrix.read_data(data)

# get probabilities each element finds exactly one match
rows_nonzero = [0] * nrow
cols_nonzero = [0] * ncol
for i in range(nrow):
    for j in range(ncol):
        if matrix.get((i, j)) > 0:
            rows_nonzero[i] += 1
            cols_nonzero[j] += 1
rows_nonzero = array(rows_nonzero)
cols_nonzero = array(cols_nonzero)

rows_optclust = array([ncol] * len(rows_nonzero)) / rows_nonzero
cols_optclust = array([nrow] * len(cols_nonzero)) / cols_nonzero

# get expected number of nonzero elements in that one cluster
# it's just 1 since we calculated # of clusters to make one cluster center with the nonzero element

rows_weighted_optclust = float(rows_optclust.dot(rows_nonzero - 1)) / sum(rows_nonzero - 1)
cols_weighted_optclust = float(cols_optclust.dot(cols_nonzero - 1)) / sum(cols_nonzero - 1)

print rows_weighted_optclust, cols_weighted_optclust
for i in range(1, 300):
    print rows_weighted_optclust / i, cols_weighted_optclust / i

fout = open("/Users/beth/Desktop/rows.txt", "w")
for i in range(len(rows_optclust)):
    print >> fout, i, rows_optclust[i]
fout.flush()
fout.close()

fout = open("/Users/beth/Desktop/cols.txt", "w")
for i in range(len(cols_optclust)):
    print >> fout, i, cols_optclust[i]
fout.flush()
fout.close()