from random import shuffle

from numpy import mean

from matrix import SparseMatrix

nrow = 3514
ncol = 1232

data = [l.split('\t') for l in open('../resources/matrix-ebc-paper-dense.tsv', 'r').readlines()]
data1 = [[e[0], e[2], e[4]] for e in data]
matrix = SparseMatrix([nrow, ncol])
matrix.read_data(data1)

# get probabilities each row in same cluster
rows = []
for i in range(nrow):
    row = []
    for j in range(ncol):
        row.append(matrix.get((i, j)))
    rows.append(tuple(row))

fout = open("/Users/beth/Desktop/test.txt", "w")
for k in range(1, nrow + 1):
    ndups = []
    for t in range(20):
        shuffle(rows)
        row_choices = rows[:k]
        other_rows = rows[k:(k + 40)]
        n_exact_matches = 0
        for row in other_rows:
            count_centers = row_choices.count(row)
            if count_centers == 1:
                n_exact_matches += 1
        ndups.append(float(n_exact_matches) / len(other_rows))
    print >> fout, k, mean(ndups)
    fout.flush()
fout.flush()
fout.close()
