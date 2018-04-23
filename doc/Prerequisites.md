# Prerequisites

## Hardware
We have tested the python components which we have written on
* Windows (Win 7 x64)
* Ubuntu 16.04 LTS.
Both platforms had >=16GB RAM. 

Upstream bioinformatics software (Kraken, Stampy, Bowtie, samtools/bcftools) have been tested on Ubuntu linux 16.04LTS.

## Software
### External tools
[Kraken v1.0](https://ccb.jhu.edu/software/kraken/), a metagenomic classifier  
[Stampy v1.0](http://www.well.ox.ac.uk/project-stampy), a mapper  
[Bowtie v2.3.4.1](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1), an alternative mapper  
[bcftools](http://www.htslib.org/doc/bcftools.html), for vcf/bcf manipulation   
[samtools](http://www.htslib.org/doc/samtools.html), for vcf/bcf manipulation  
[python v3.5.2](https://www.python.org/downloads/)

Notes:  
* a script to obtain these dependencies is [available](../pipeline/get_software.sh) (see also below).
* the test data set does not require any tools, except python, to be on the PATH
* The software is not expected to work with python 2.7.
* Although the paper used R from some analyses, in these have all been ported to Python in the software released here.

### Python packages
The following modules which are not part of the python3 standard library are required

[numpy](http://www.numpy.org/) (tested with 1.14.0)  
[scipy](https://www.scipy.org/) (tested with 1.0.0)  
[pandas](https://pandas.pydata.org/) (tested with 0.22.0)  
[tables](https://www.pytables.org/) (tested with 3.1)  
[Biopython](http://biopython.org/) (tested with 1.7.0)

If they are not present, they can be installed using commands similar to
```
sudo pip3 install tables
```

*for fitting regression models*  
[statsmodels](http://www.statsmodels.org/stable/index.html) (tested with 0.8)

*for visualisation of output*  
[bokeh](https://bokeh.pydata.org/en/latest/) (tested with 0.12.14)

On Windows systems, we used the numpy+mkl windows binaries [available here](https://www.lfd.uci.edu/~gohlke/pythonlibs/).
On both platforms, we installed all other dependencies using pip3.

## Getting software for an end-to-end demonstration of the process
This can be installed automatically on linux systems.  Note:
* The paths given for software downloads may change as third party software versions update.  Please see the software URLs above ('software') if links break.
* The below script will also download a 4GB minikraken database.  Please edit the script if you do not require this.

```
cd pipeline/software
chmod +x get_software
./get_software.sh
```

Output is written into the /pipeline/software directory.

## Test data
Some test data is supplied when cloning the github repository, but 
* vcf files used for unit testing
* fastq.gz files used as input to the end-to-end demonstration
  are too large to be installed from github.

This data is being transferred to a permanent repository  [here](https://ora.ox.ac.uk/objects/uuid:88d93ec1-2757-4e46-83f4-dbdfc01b5343)
but is currently available [here](https://www.dropbox.com/sh/rblwq86qd6cnjfb/AAB35TApL8KLfJINReTPCSsya).
You can obtain test data as below:

### test data for unit testing [1.2GB]
```
cd testdata
DOWNLOADURL=https://www.dropbox.com/s/277ops3cg2ngu2k/testdata.tar
wget $DOWNLOADURL
tar -xf testdata.tar
```

### input data comprising 250 fastq.gz files for an end-to-end demonstration [65G]
```
cd pipeline/testdata
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_0.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_1.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_2.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_3.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_4.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_5.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_6.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_7.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_8.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_9.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_a.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_b.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_c.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_d.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_e.tar
wget https://www.dropbox.com/s/vhp7rv6a3hpqz2j/fq_f.tar

tar -xf *.tar

```
