from collections import defaultdict
import os

matrix2d = defaultdict(int)
matrix3d = defaultdict(int)

file_count = 0
for root, dirs, files in os.walk("/rd02/data/blpercha/ebc-results/"):
    for f in files:
        print file_count, f
        file_count += 1
        file_open = open(os.path.join(root, f), "r")
        if "3d" in f:
            for line in file_open:
                sl = line.split("\t")
                p1 = "(" + ",".join(sl[2].split(" ")) + ")"
                p2 = "(" + ",".join(sl[3].split(" ")) + ")"
                count = int(sl[4])
                ids = p1 + "\t" + "-" + "\t" + p2 + "\t" + "-"
                matrix3d[ids] += count
        else:
            for line in file_open:
                sl = line.split("\t")
                p1 = sl[2]
                p2 = sl[3]
                count = int(sl[4])
                ids = p1 + "\t" + "-" + "\t" + p2 + "\t" + "-"
                matrix2d[ids] += count

output_2d = open("/home/blpercha/ebc-results-combined-2d.txt", "w")
for k in matrix2d:
    output_2d.write(k + "\t" + str(matrix2d[k]) + "\n")
output_2d.flush()
output_2d.close()

output_3d = open("/home/blpercha/ebc-results-combined-3d.txt", "w")
for k in matrix3d:
    output_3d.write(k + "\t" + str(matrix3d[k]) + "\n")
output_3d.flush()
output_3d.close()