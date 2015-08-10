f = open("resources/matrix-ebc-paper-sparse.tsv", "r")
f3d = open("resources/matrix-ebc-paper-sparse-3d.tsv", "w")

for line in f:
    sl = line.split("\t")
    col1_names = sl[0].strip("()").split(",")
    col2_name = sl[2]
    val = float(sl[4].strip())
    print >> f3d, col1_names[0] + "\t" + col1_names[1] + "\t" + col2_name + "\t" + str(val)

f.close()
f3d.flush()
f3d.close()