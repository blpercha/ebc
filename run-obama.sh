#!/usr/bin/env bash

 for i in {1..200}
    do bsub -M 2 -R "rusage[mem=2]" \
        -o /home/blpercha/output-logs/${i}.out -e /home/blpercha/output-logs/${i}.err \
        python /home/blpercha/ebc/clusters.py /home/blpercha/ebc/resources/matrix-ebc-paper-dense.tsv \
        0,2,4 3000 5 10 /home/blpercha/ebc-results/matrix-ebc-paper-dense-clustopt-${i}.txt
 done
