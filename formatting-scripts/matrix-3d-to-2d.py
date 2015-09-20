import sys

f = open(sys.argv[1], "r")
f3d = open(sys.argv[2], "w")

for line in f:
    sl = line.split("\t")
    col1_name1 = sl[0]
    col1_name2 = sl[1]
    col2_name = sl[2]
    val = float(sl[3].strip())
    print >> f3d, "(" + ",".join([col1_name1, col1_name2]) + ")" + "\t" + col2_name + "\t" + str(val)

f.close()
f3d.flush()
f3d.close()
