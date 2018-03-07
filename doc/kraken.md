# Running Kraken

[Kraken](https://ccb.jhu.edu/software/kraken/) is a metagenomic classifier which assigns reads to taxa based on a dictionary.
For the publication, we used a custom database constructed from bacterial reference sequences in Genbank, plus human sequences.
Instructions on how to build this are in Methods.

For this demonstration, we illustrate the process using one of the [MiniKraken databases](https://ccb.jhu.edu/software/kraken/dl/minikraken_20171101_4GB_dustmasked.tgz) released by Kraken's developers.
This has similar specificity, but about 25% lower sensitivity, for genus-level read assignment.

Example usage is as follows:

```
kraken --threads 1 --preload --db /mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25 --fasta-input /dev/fd/0
kraken-filter.o --threads 1 --db /mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25 --threshold 0.05
kraken-report.o --db /mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25
```
