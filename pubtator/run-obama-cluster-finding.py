import os
import sys

matrix_file = sys.argv[1]
output_file = sys.argv[2]
memory = int(sys.argv[3])

for i in xrange(5, 100, 10):
    for j in xrange(105, 200, 10):
        output_log_file = str(i) + "-" + str(j)
        os.system("bsub " +
                  "-o /home/blpercha/output-logs/" + output_log_file + ".out " +
                  "-e /home/blpercha/output-logs/" + output_log_file + ".err " +
                  "-M %d -R \"rusage[mem=%d]\" " % (memory, memory) +
                  "python /home/blpercha/ebc/clusters-single.py " +
                  "%s " % matrix_file +
                  "0,1,2 " +
                  "%d,%d " % (i, j) + "20 %s 1e-10 2" % output_file)
