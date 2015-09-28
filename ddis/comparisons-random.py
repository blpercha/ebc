from collections import defaultdict
import csv

# get known ddis
from math import log
import os
from random import shuffle
from scipy.stats import hypergeom
from sklearn import metrics

known_ddis = set()
f2 = csv.reader(open("/Users/beth/Documents/phd/ebc-ddis/pgkb-druginteractions-clean-20110414.csv", "r"),
                delimiter=',', quotechar='\"')

for line in f2:
    d1 = line[1]
    d2 = line[3]
    known_ddis.add((d1, d2))
    known_ddis.add((d2, d1))

print "Got known ddis"


# collect all drug pairs in matrix
all_drug_pairs = set()
clusters = defaultdict(list)
cluster_counts = defaultdict(lambda: defaultdict(int))

f3d = open("/Users/beth/Desktop/ebc-results-combined-ddis-iterations-2.tsv", "r")

for line in f3d:
    sl = line.split("\t")
    drug_pair = tuple(sl[0].split(","))
    all_drug_pairs.add(drug_pair)
    clusters[drug_pair] = [int(e) for e in sl[1:]]
    for t, c in enumerate(clusters[drug_pair]):
        cluster_counts[t][c] += 1

f3d.close()

known_ddis_matrix = list(all_drug_pairs & known_ddis)
non_ddis_matrix = list(all_drug_pairs - known_ddis)

print "Number of known ddis in matrix: ", len(known_ddis_matrix)
print "Number of additional rows in matrix: ", len(non_ddis_matrix)

# make training and test sets and run comparisons
test_set_size = 100

output_folder = "/Users/beth/Desktop/ddi-results"

M = len(clusters)  # total number of rows in matrix i.e. number of drug pairs
for S in [1, 2, 3, 4, 5, 10, 25, 50, 75, 100]:
    fout = open(os.path.join(output_folder, str(S) + ".txt"), "w")
    for j in range(100):
        shuffle(known_ddis_matrix)
        shuffle(non_ddis_matrix)

        # figure out which clusters the seed members are in
        seed = known_ddis_matrix[:S]
        cluster_seed_counts = defaultdict(lambda : defaultdict(int))
        for s in seed:
            for t in range(len(clusters[s])):
                cluster_seed_counts[t][clusters[s][t]] += 1

        # calculate cluster scores based on the hypergeometric distribution
        cluster_scores = defaultdict(lambda: defaultdict())
        for t in cluster_counts:
            for c in cluster_counts[t]:
                N = cluster_counts[t][c]
                x = cluster_seed_counts[t][c]
                # print x, M, S, N, 1.0 - hypergeom.cdf(x, M, S, N)
                cluster_scores[t][c] = 1.0 - hypergeom.cdf(x, M, S, N)

        test_positive = known_ddis_matrix[S:(S + test_set_size / 2)]
        test_negative = non_ddis_matrix[:(test_set_size / 2)]

        # compare test data with seed set
        scores = []
        true_labels = []
        for x in test_positive:
            score = 0.0
            for t, c in enumerate(clusters[x]):
                # print cluster_scores[t][c]
                score += cluster_scores[t][c]
            scores.append(score)
            true_labels.append(1)
        for x in test_negative:
            score = 0.0
            for t, c in enumerate(clusters[x]):
                score += cluster_scores[t][c]
            scores.append(score)
            true_labels.append(0)

        fpr, tpr, _ = metrics.roc_curve(true_labels, scores)
        auc = metrics.auc(fpr, tpr)
        fout.write(str(auc) + "\n")
        fout.flush()
    fout.close()
