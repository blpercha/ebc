import os

k_0 = 35
k_1 = 15
k_2 = 15
n_runs = 50
n_simultaneous = 40
n_iterations = 2

for t in range(n_simultaneous):
    output_log_file = "repeat-runs-ddis"
    os.system("bsub -M 1 -R \"rusage[mem=1]\" " +
              "-o /home/blpercha/output-logs/" + output_log_file + ".out " +
              "-e /home/blpercha/output-logs/" + output_log_file + ".err " +
              "python /home/blpercha/ebc/clusters-repeat.py " +
              "/home/blpercha/matrix-ebc-paper-ddi-format.tsv " +
              "0,3,4,6 " +
              "%d,%d,%d " % (k_0, k_1, k_2) +
              "%d " % n_runs +
              "/rd02/data/blpercha/ebc-results-ddis/repeat-%d.txt " % t +
              "1e-10 " +
              "%d 0" % n_iterations)
