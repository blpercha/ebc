import sys

f = open(sys.argv[1], "r")
f3d = open(sys.argv[2], "w")

for line in f:
    sl = line.split("\t")
    names = sl[:(len(sl) - 1)]
    print >> f3d, "\t".join(names) + "\t1.0"

f.close()
f3d.flush()
f3d.close()
