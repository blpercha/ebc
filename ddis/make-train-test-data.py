from collections import defaultdict
import csv

# get known ddis
import os
from random import shuffle
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

f3d = open("/Users/beth/Desktop/ebc-results-combined-ddis-iterations-2.tsv", "r")

for line in f3d:
    sl = line.split("\t")
    drug_pair = tuple(sl[0].split(","))
    all_drug_pairs.add(drug_pair)
    clusters[drug_pair] = [int(e) for e in sl[1:]]

f3d.close()

print "Got all drug pairs"


# make training and test sets and run comparisons

def compare(v1, v2):
    same = 0
    for i in range(len(v1)):
        if v1[i] == v2[i]:
            same += 1
    return same


known_ddis_matrix = list(all_drug_pairs & known_ddis)
non_ddis_matrix = list(all_drug_pairs - known_ddis)

test_set_size = 100

output_folder = "/Users/beth/Desktop/ddi-results"

for i in [1, 2, 3, 4, 5, 10, 25, 50, 75, 100]:
    fout = open(os.path.join(output_folder, str(i) + ".txt"), "w")
    for j in range(100):
        shuffle(known_ddis_matrix)
        shuffle(non_ddis_matrix)

        seed = known_ddis_matrix[:i]
        test_positive = known_ddis_matrix[i:(i + test_set_size / 2)]
        test_negative = non_ddis_matrix[:(test_set_size / 2)]

        # compare test data with positive training data
        scores = []
        true_labels = []
        for x in test_positive:
            sum_seed_scores = 0.0
            for t in seed:
                sum_seed_scores += float(compare(clusters[x], clusters[t]))
            scores.append(sum_seed_scores)
            true_labels.append(1)
        for x in test_negative:
            sum_seed_scores = 0.0
            for t in seed:
                sum_seed_scores += float(compare(clusters[x], clusters[t]))
            scores.append(sum_seed_scores)
            true_labels.append(0)

        fpr, tpr, _ = metrics.roc_curve(true_labels, scores)
        auc = metrics.auc(fpr, tpr)
        fout.write(str(auc) + "\n")
        fout.flush()
    fout.close()
