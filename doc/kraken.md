# Running Kraken

[Kraken](https://ccb.jhu.edu/software/kraken/) is a metagenomic classifier which assigns reads to taxa based on a dictionary.  
For the publication, we used a custom database constructed from bacterial RefSeq in Genbank, plus a  human genomic build.  
Instructions on how to build this database are in Methods.

For this demonstration, we illustrate the process using one of the [MiniKraken databases](https://ccb.jhu.edu/software/kraken/dl/minikraken_20171101_4GB_dustmasked.tgz) released by Kraken's developers.  
Their data indicates that compared with the full database, this mini-database has similar specificity, but about 25% lower sensitivity, for genus-level read assignment.

Example usage is as follows:

```
kraken --threads 1 --preload --db /mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25 --fasta-input /dev/fd/0
kraken-filter.o --threads 1 --db /mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25 --threshold 0.05
kraken-report.o --db /mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25
```

The KrakenReportReader class will summarise data from Kraken/Braken reports.
Please see the documentation in the KrakenReportReader module for more information.

```python
from KrakenReportReader import KrakenReportReader

# example usage
krr = KrakenReportReader(genus_of_interest = 'Mycobacterium')
inputfile = os.path.join('..','testdata','test1.kraken_report')
result = krr.simplify(inputfile= inputfile, guid = 'test1')


```

### Command line usage
See above for detail of how to use kraken.  
The *KrakenReportReader* class is used by the *AdaptiveMasking* class.  You do not need to use it directly. See [here](model_maf.md).
