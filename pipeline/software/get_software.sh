# get software
echo Getting kraken
wget https://ccb.jhu.edu/software/kraken/dl/kraken-1.0.tgz
tar -xvf kraken-1.0.tgz

echo Getting Samtools
wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
tar -xvf samtools-1.7.tar.bz2

echo Getting BcfTools
wget https://github.com/samtools/bcftools/releases/download/1.7/bcftools-1.7.tar.bz2
tar -xvf bcftools-1.7.tar.bz2

echo Getting GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.0.2.1/gatk-4.0.2.1.zip
unzip gatk-4.0.2.1.zip

echo Getting Bowtie
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-linux-x86_64.zip
unzip bowtie-1.2.2-linux-x86_64.zip

echo Getting Stampy
wget http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz
tar -xvf Stampy-latest.tgz

echo Getting MiniKraken dusted database
wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171101_4GB_dustmasked.tgz
tar -xvf minikraken_20171101_4GB_dustmasked.tgz

echo Cleanup
rm *.zip
rm *.tgz
rm *.gz2

echo Compiling Stampy
cd stampy-1.0.32
make
echo Compiling Bcftools
cd ..
cd bcftools-1.7
make
echo Compliling Samtools
cd  ..
cd samtools-1.7
make
