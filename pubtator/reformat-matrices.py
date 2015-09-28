import os

for root, dirs, files in os.walk("/Users/beth/Desktop/pubtator-ebc-matrices-dependency-paths"):
    for f in files:
        if ".txt" not in f:
            continue

        x = open(os.path.join(root, f), "r")
        y = open(os.path.join(root, f + "2"), "w")
        for line in x:
            sl = line.split("\t")
            col1_name1 = sl[0]
            col1_name2 = sl[1]
            if col1_name1 == col1_name2:
                continue
            col2_name = sl[2]
            val = float(sl[3].strip())
            y.write("(" + ",".join([col1_name1, col1_name2]) + ")" + "\t" + col2_name + "\t" + str(1.0) + "\n")

