import os

for i in xrange(10, 150, 10):
    for j in xrange(10, 150, 10):
        for k in xrange(10, 150, 10):
            output_log_file = str(i) + "-" + str(j) + "-" + str(k)
            os.system("bsub -M 1 -R \"rusage[mem=1]\" " +
                      "-o /home/blpercha/output-logs/" + output_log_file + ".out " +
                      "-e /home/blpercha/output-logs/" + output_log_file + ".err " +
                      "python /home/blpercha/ebc/clusters-single.py " +
                      "/home/blpercha/ebc/resources/matrix-ebc-paper-dense-3d.tsv " +
                      "0,1,2,3 " +
                      "%d,%d,%d " % (i, j, k) + "20 /home/blpercha/ebc-results-3d.txt")
