import os

for i in xrange(80, 100, 2):
    for j in xrange(100, 120, 2):
        for k in xrange(2, 16, 2):
            output_log_file = str(i) + "-" + str(j) + "-" + str(k)
            os.system("bsub -M 1 -R \"rusage[mem=1]\" " +
                      "-o /home/blpercha/output-logs/" + output_log_file + ".out " +
                      "-e /home/blpercha/output-logs/" + output_log_file + ".err " +
                      "python /home/blpercha/ebc/clusters-single.py " +
                      "/home/blpercha/ebc/resources/matrix-ebc-paper-dense-3d.tsv " +
                      "0,1,2,3 " +
                      "%d,%d,%d " % (i, j, k) + "50 /home/blpercha/ebc-results-3d-gridsize-2.txt 1e-10 30")
