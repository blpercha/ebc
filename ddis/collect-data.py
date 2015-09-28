from collections import defaultdict
import os

final_data = defaultdict(list)

count = 0
for root, dirs, files in os.walk("/rd02/data/blpercha/ebc-results-ddis/"):
    for fil in files:
        print count, fil
        f = open(os.path.join(root, fil), "r")
        for line in f:
            sl = line.split("\t")
            id_pair = sl[0]
            names = sl[1]
            clusters = [int(e) for e in sl[2:]]
            for c in clusters:
                final_data[names].append(c)
        count += 1

# write final data file

fout = open("/home/blpercha/ebc-results-combined-ddis-iterations-2.tsv", "w")

for names in final_data:
    print >> fout, names + "\t" + "\t".join([str(e) for e in final_data[names]])

fout.flush()
fout.close()