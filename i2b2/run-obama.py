import os

for i in xrange(5, 100, 10):
    for j in xrange(105, 200, 10):
            output_log_file = str(i) + "-" + str(j)
            os.system("bsub " +
                      "-o /home/blpercha/output-logs/" + output_log_file + ".out " +
                      "-e /home/blpercha/output-logs/" + output_log_file + ".err " +
                      "python /home/blpercha/ebc/clusters-single.py " +
                      "/home/blpercha/matrix-i2b2-words.txt " +
                      "0,1,2 " +
                      "%d,%d " % (i, j) + "20 /home/blpercha/ebc-results-i2b2-gridsize-10.txt 1e-10 2")
