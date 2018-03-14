# Running Kraken

## Background
[Kraken](https://ccb.jhu.edu/software/kraken/) is a metagenomic classifier which assigns reads to taxa based on a dictionary.  
For the publication, we used a custom database constructed from bacterial RefSeq in Genbank, plus a  human genomic build.  
Instructions on how to build this database are in Methods.

For this demonstration, we illustrate the process using one of the [MiniKraken databases](https://ccb.jhu.edu/software/kraken/dl/minikraken_20171101_4GB_dustmasked.tgz) released by Kraken's developers.  
Their data indicates that compared with the full database, this mini-database has similar specificity, but about 25% lower sensitivity, for genus-level read assignment.

Example usage is as follows:

```
kraken --threads 1 \
--preload --db /mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25 \
--fastq-input --paired testdata/in.r1.fq testdata/in.r2.fq |
kraken-filter --db \
/mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25 --threshold 0.05 > testdata/in.kraken.txt

kraken-report --db /mnt/microbio/pipeline/KRAKEN/kraken_full_18012016_added_25 \
testdata/in.kraken.txt > testdata/in.kraken.report.txt
```

## Pipeline
For a further example, please see *pipeline/testdata/run_kraken.sh*.  This script runs Kraken against fifty fastq files,
obtained from a clinical *M. tuberculosis* sequencing operation.  

*1 Obtain the test data and software*   
The software required and test data can be obtained as [described](Prerequisites.md)

*2 Run mapping*    
A shell script is provided to do this.  

```
cd pipeline/testdata
chmod +x run_kraken.sh

# recommend running with nohup, as it will take 15-30 minutes depending on the architecture used.
# note: the shell script checks for output files; it will not rerun analyses
nohup ./run_kraken.sh > krak.out 2> krak.err &
watch krak.out

```


## Summarising the output
The KrakenReportReader class will summarise data from Kraken/Braken reports.
Please see the documentation in the KrakenReportReader module for more information.

```python
from KrakenReportReader import KrakenReportReader

# example usage
krr = KrakenReportReader(genus_of_interest = 'Mycobacterium')
inputfile = os.path.join('..','testdata','test1.kraken_report')
result = krr.simplify(inputfile= inputfile, guid = 'test1')


```

The *KrakenReportReader* class is used by the *AdaptiveMasking* class.  You do not need to use it directly. See [here](model_maf.md).
