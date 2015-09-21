import os
import sys
from collections import defaultdict


# args: 1 - directory with results files
#       2 - output combined file

fout = open(sys.argv[2], "w")

full_counts = defaultdict(int)

count = 1
for root, dirs, files in os.walk(sys.argv[1]):
    for fil in files:
        if ".txt" not in fil:
            continue
        print(count, fil)
        count += 1
        f = open(os.path.join(root, fil), "r")
        for line in f:
            split_line = line.split("\t")
            ids = "\t".join([split_line[2], split_line[0], split_line[3], split_line[1]])
            value = int(split_line[4])
            full_counts[ids] += value
        f.close()

for k in full_counts:
    fout.write("%s\t%d\n" % (k, full_counts[k]))

fout.flush()
fout.close()
