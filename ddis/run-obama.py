import os

for i in xrange(5, 50, 10):
    for j in xrange(5, 50, 10):
        for k in xrange(5, 50, 10):
            output_log_file = str(i) + "-" + str(j) + "-" + str(k)
            os.system("bsub " +
                      "-o /home/blpercha/output-logs/" + output_log_file + ".out " +
                      "-e /home/blpercha/output-logs/" + output_log_file + ".err " +
                      "python /home/blpercha/ebc/clusters-single.py " +
                      "/home/blpercha/matrix-ebc-paper-ddi-format.tsv " +
                      "0,3,4,6 " +
                      "%d,%d,%d " % (i, j, k) + "20 /home/blpercha/ebc-results-ddis-gridsize-10.txt 1e-10 100")
