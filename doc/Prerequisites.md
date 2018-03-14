# Prerequisites

## Hardware
We have tested the python components which we have written on
* Windows (Win 7 x64)
* Ubuntu 16.04 LTS.
Both platforms had >=16GB RAM. 

Upstream bioinformatics software (Kraken, Stampy, Bowtie, samtools/bcftools) have only been tested on linux.

## Software
### External tools
[Kraken v1.0](https://ccb.jhu.edu/software/kraken/), a metagenomic classifier  
[Stampy v1.0](http://www.well.ox.ac.uk/project-stampy), a mapper  
[Bowtie v2.3.4.1](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1), an alternative mapper  
[bcftools](http://www.htslib.org/doc/bcftools.html), for vcf/bcf manipulation   
[samtools](http://www.htslib.org/doc/samtools.html), for vcf/bcf manipulation
bedtools (for bam to fastq conversion, if mapping using test data).  The test scripts assume this is on the path.
[python v3.5.2](https://www.python.org/downloads/)      


The software is not expected to work with python 2.7.
Although the paper used R from some analyses, in these have all been ported to Python in the software released here.

### Python packages
The following modules which are not part of the python3 standard library are required

[numpy](http://www.numpy.org/) (tested with 1.14.0)  
[scipy](https://www.scipy.org/) (tested with 1.0.0)  
[pandas](https://pandas.pydata.org/) (tested with 0.22.0)  
[tables](https://www.pytables.org/) (tested with 3.1)  
[Biopython](http://biopython.org/) (tested with 1.7.0)

*for fitting regression models*  
[statsmodels](http://www.statsmodels.org/stable/index.html) (tested with 0.8)

*for visualisation of output*  
[bokeh](https://bokeh.pydata.org/en/latest/) (tested with 0.12.14)


On Windows systems, we used the numpy+mkl windows binaries [available here](https://www.lfd.uci.edu/~gohlke/pythonlibs/).
On both platforms, we installed all other dependencies using pip3.

## Getting software
This can be installed automatically on linux systems.  The below script will also download a minikraken database.
```
cd pipeline/software
chmod +x get_software
./get_software.sh
```

See the /pipeline/software directory
## Test data
Some test data is supplied when cloning the github repository, but some test files (including vcf files) are too large to be installed from github.
You can obtain this data as below:
```
cd testdata
wget ****
gunzip *

```
You can also do so by running the get_test_data.py script in the src/ directory.

```
cd src
python3 get_test_data.py
```
