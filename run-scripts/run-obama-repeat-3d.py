import os

k_0 = 98
k_1 = 112
k_2 = 2
n_runs = 50
n_simultaneous = 40

for t in range(n_simultaneous):
    output_log_file = str(t) + ".txt"
    os.system("bsub -M 1 -R \"rusage[mem=1]\" " +
              "-o /home/blpercha/output-logs/" + output_log_file + ".out " +
              "-e /home/blpercha/output-logs/" + output_log_file + ".err " +
              "python /home/blpercha/ebc/clusters-repeat.py " +
              "/home/blpercha/ebc/resources/matrix-ebc-paper-dense-3d.tsv " +
              "0,1,2,3 " +
              "%d,%d,%d " % (k_0, k_1, k_2) +
              "%d " % n_runs +
              "/home/blpercha/ebc-results/repeat-3d-%d.txt " % t +
              "1e-10 " +
              "100 0,1")
