# Prerequisites

## Hardware
We have tested the R and python components which we have written on
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
[python v3.5.2](https://www.python.org/downloads/)      
[R v3.3](https://cran.r-project.org/)    
 
No testing has been done with earlier versions.  It is unlikely the software will work with python 2.7.

### Python packages
The following modules which are not part of the python3 standard library are required

numpy  
scipy  
pandas  
tables  
Biopython  

On Windows systems, we used numpy+mkl windows binaries [available here](https://www.lfd.uci.edu/~gohlke/pythonlibs/)
and installed other dependencies using pip3.

## R modules
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
