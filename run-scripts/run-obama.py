import os
import sys

matrix_file = sys.argv[1]
output_file = sys.argv[2]
memory = int(sys.argv[3])

for i in xrange(10, 200, 10):
    for j in xrange(10, 200, 10):
        output_log_file = str(i) + "-" + str(j)
        os.system("bsub -M %d -R \"rusage[mem=%d]\" " % (memory, memory) +
                  "-o /home/blpercha/output-logs/" + output_log_file + ".out " +
                  "-e /home/blpercha/output-logs/" + output_log_file + ".err " +
                  "python /home/blpercha/ebc/clusters-single.py " +
                  "%s " % matrix_file +
                  "0,1,2 " +
                  "%d,%d " % (i, j) + "10 %s 1e-10 50" % output_file)
