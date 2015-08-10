#!/usr/bin/env bash

 for i in {1..400}
    do bsub -M 1 -R "rusage[mem=1]" \
        -o /home/blpercha/output-logs/${i}.out -e /home/blpercha/output-logs/${i}.err \
        python /home/blpercha/ebc/clusters.py /home/blpercha/ebc/resources/matrix-ebc-paper-dense.tsv \
        0,2,4 200 10 5 /home/blpercha/ebc-results/matrix-ebc-paper-dense-clustopt-${i}.txt 10.0
 done
