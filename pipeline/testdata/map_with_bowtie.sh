#!/bin/sh
# start in testdata;
cd ../testdata

# construct paths to executables relative to the pipeline directory
BOWTIE=software/bowtie-1.2.2-linux-x86_64/
STAMPY=software/stampy-1.0.32/
SAMTOOLS=software/samtools-1.7/#
# define reference fasta file for readability
REFFASTA=../../../testdata/NC_000962.3.fasta
for dir in *
do
    GUID=${dir%*/}
    if [ -d $GUID ]; then
        echo Analysing: ${GUID}     
        cd ${GUID}
                bedtools bamtofastq -i ${GUID}.bam -fq ${GUID}.fq
        gzip -f ${GUID}.fq
        ../../${BOWTIE}bowtie ../../mapto/tb_bowtie ${GUID}.fq.gz bowtie_output.sam -S
        ../../${SAMTOOLS}samtools view -t ../../../testdata/samtools_tab_file.txt  -o bowtie_output.bam -b -S bowtie_output.sam
        rm -f bowtie_output.sam
        ../../${SAMTOOLS}samtools sort bowtie_output.bam -t sort_tmpfile -o sorted_bowtie_output.bam -m 8G
        rm -f bowtie_output.bam
        ../../${SAMTOOLS}samtools index sorted_bowtie_output.bam
        # can compute depth from INFO/AD tag
        ../../${SAMTOOLS}samtools mpileup -A -B -f ${REFFASTA} -Q25 -q30 -o40 -e20 -h100 -m2 -F0.002  -t INFO/AD -s -v  sorted_bowtie_output.bam -o ${GUID}.vcf.gz
        rm -f sorted_bowtie_output.bam
        rm -f sorted_bowtie_output.bam.bai
        cd ..
    fi
done
exit 0



