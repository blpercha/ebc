f = open("../resources/matrix-ebc-paper-dense.tsv", "r")
f4d = open("../resources/matrix-ebc-paper-dense-4d.tsv", "w")

for line in f:
    sl = line.split("\t")
    col1_names_joined = sl[0]
    col1_names = sl[0].strip("()").split(",")
    col2_name = sl[2]
    val = float(sl[4].strip())
    print >> f4d, col1_names_joined + "\t" + col1_names[0] + "\t" + col1_names[1] + "\t" + \
                  col2_name + "\t" + str(val)

f.close()
f4d.flush()
f4d.close()
