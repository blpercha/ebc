import os

for root, dirs, files in os.walk("/Users/beth/Desktop/pubtator-ebc-matrices"):
    for f in files:
        if ".txt" not in f:
            continue
        x = open(os.path.join(root, f), "r")
        y = open(os.path.join(root, f + "2"), "w")
        for line in x:
            sl = line.split("\t")
            if sl[0] == sl[1]:
                continue
            y.write(line)
