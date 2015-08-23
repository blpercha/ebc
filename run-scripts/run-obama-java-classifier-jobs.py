import os
import sys

ebc_same_cluster_file = sys.argv[1]
raw_matrix = sys.argv[2]
training_folder_stem = sys.argv[3]
test_folder_stem = sys.argv[4]
output_folder_stem = sys.argv[5]
memory = int(sys.argv[6])

num_datapoints = [1, 2, 3, 4, 5, 10, 25, 50, 75, 100]

os.system("module load java8")

for ndat in num_datapoints:
    training_folder = os.path.join(training_folder_stem, "train-%d" % ndat)
    test_folder = os.path.join(test_folder_stem, "test-%d" % ndat)

    # submit ebcranksum jobs
    output_folder = os.path.join(output_folder_stem, "results-%d-ebcranksum" % ndat)
    os.system("rm -r %s" % output_folder)
    os.system("mkdir %s" % output_folder)
    for root, dirs, files in os.walk(training_folder):
        for training_file in files:
            if ".tsv" not in training_file:
                continue
            training_file_full = os.path.join(root, training_file)
            file_base = os.path.splitext(training_file)[0]
            test_file = os.path.join(test_folder, training_file)
            output_file = os.path.join(output_folder, file_base + ".tsv")
            os.system(
                "bsub -M %d -R \"rusage[mem=%d]\" -o output-logs/%s.out -e output-logs/%s.err java -mx%dg -cp /home/blpercha/reflective-dirt.jar classifier.RunClassifierMain -t ebcranksum -r %s %s %s %s" % (
                    memory, memory, file_base, file_base, memory, ebc_same_cluster_file, training_file_full, test_file,
                    output_file))
