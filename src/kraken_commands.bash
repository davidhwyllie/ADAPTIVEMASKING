#!/bin/bash

# example commands for guuid 745392e5-0a94-4f0d-b60a-7888312f1729
# converted from bam to fastq usibng bedtools
# example: bamToFastq -i in.bam -fq in.r1.fq -fq2 in.r2.fq
kraken --threads 1 \
--preload --db /mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25 \
--fastq-input --paired testdata/in.r1.fq testdata/in.r2.fq |
kraken-filter --db \
/mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25 --threshold 0.05 > testdata/in.kraken.txt

kraken-report --db /mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25 \
testdata/in.kraken.txt > testdata/in.kraken.report.txt
