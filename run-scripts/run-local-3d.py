import os

for i in xrange(1, 6, 1):
    for j in xrange(1, 6, 1):
        for k in xrange(1, 6, 1):
            output_log_file = str(i) + "-" + str(j) + "-" + str(k)
            os.system("python /Users/beth/Documents/phd/ebc/clusters-single.py " +
                      "/Users/beth/Documents/phd/ebc/resources/matrix-itcc-3d-3clust.tsv " +
                      "0,1,2,3 " +
                      "%d,%d,%d " % (i, j, k) + "50 /Users/beth/Desktop/ebc-results-3d-gridsize-2.txt 1e-10 100")
