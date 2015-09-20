from collections import defaultdict
import os


def get_condition_labels(f):
    labels = defaultdict(list)
    for line in f:
        sl = line.split("\t")
        if sl[2].strip() == 'Y':
            condition = sl[1]
            note = sl[0]
            labels[condition].append(note)
    return labels


labels_intuitive = get_condition_labels(open(
    "/Users/beth/Documents/phd/ebc-i2b2/i2b2-2008-notes/all_intuitive.txt", "r"))

labels_textual = get_condition_labels(open(
    "/Users/beth/Documents/phd/ebc-i2b2/i2b2-2008-notes/all_textual.txt", "r"))

# get words from notes and write as matrix
fout = open("/Users/beth/Documents/phd/ebc-i2b2/matrix-i2b2-words.txt", "w")
for root, dirs, files in os.walk("/Users/beth/Documents/phd/ebc-i2b2/i2b2-2008-notes"):
    for fil in files:
        if "note" not in fil:
            continue
        fil_open = open(os.path.join(root, fil), "r")
        words = set(fil_open.read().split())
        for w in words:
            print >> fout, fil + "\t" + w + "\t1"
fout.flush()
fout.close()
