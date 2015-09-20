import os

for i in xrange(70, 111, 5):
    for j in xrange(2, 15, 2):
        for k in xrange(5, 25, 2):
            for l in xrange(2, 15, 2):
                output_log_file = str(i) + "-" + str(j) + "-" + str(k) + "-" + str(l)
                os.system("bsub -M 1 -R \"rusage[mem=1]\" " +
                          "-o /home/blpercha/output-logs/" + output_log_file + ".out " +
                          "-e /home/blpercha/output-logs/" + output_log_file + ".err " +
                          "python /home/blpercha/ebc/clusters-single.py " +
                          "/home/blpercha/ebc/resources/matrix-ebc-paper-dense-4d.tsv " +
                          "0,1,2,3,4 " +
                          "%d,%d,%d,%d " % (i, j, k, l) + "50 /home/blpercha/ebc-results-4d-5.txt 1e-10 100")
