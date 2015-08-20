import os

for i in xrange(10, 200, 10):
    for j in xrange(10, 200, 10):
        output_log_file = str(i) + "-" + str(j)
        os.system("bsub -M 1 -R \"rusage[mem=1]\" " +
                  "-o /home/blpercha/output-logs/" + output_log_file + ".out " +
                  "-e /home/blpercha/output-logs/" + output_log_file + ".err " +
                  "python /home/blpercha/ebc/clusters-single.py " +
                  "/home/blpercha/ebc/resources/matrix-ebc-paper-dense.tsv " +
                  "0,2,4 " +
                  "%d,%d " % (i, j) + "50 /home/blpercha/ebc-results.txt 1e-10")
