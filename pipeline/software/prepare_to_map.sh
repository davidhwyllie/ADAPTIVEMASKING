# define software
# construct paths to executables relative to the pipeline directory
BOWTIE=software/bowtie-1.2.2-linux-x86_64/
STAMPY=software/stampy-1.0.32/

# create bowtie index in directory mapto
${BOWTIE}bowtie-build  ../testdata/NC_000962.3.fasta mapto/tb_bowtie
${SAMTOOLS}samtools faidx ../testdata/NC_000962.3.fasta
cd mapto
../${STAMPY}stampy.py  --species=M_tuberculosis --assembly=NC_000962.3 -G tb_stampy ../../testdata/NC_000962.3.fasta
../${STAMPY}stampy.py  -g tb_stampy -H tb_stampy
  
  