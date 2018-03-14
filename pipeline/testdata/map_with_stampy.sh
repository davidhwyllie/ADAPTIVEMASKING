!/bin/sh
# start in testdata;
cd ../testdata

# construct paths to executables relative to the pipeline directory
STAMPY=software/stampy-1.0.32/
SAMTOOLS=software/samtools-1.7/

# define reference fasta file for readability
REFFASTA=../../../testdata/NC_000962.3.fasta


for dir in *
do
    GUID=${dir%*/}
    echo Analysing: ${GUID}
    if [ -d $GUID ]; then
        cd ${GUID}
        TARGETFILE=${GUID}.fq.gz
	if [ -f $TARGETFILE ]; then
		echo 'fastq.gz read file found'
	else
		echo 'fastq.gz read file not found; generating from Bam file if present'
                bedtools bamtofastq -i ${GUID}.bam -fq ${GUID}.fq
                gzip -f ${GUID}.fq
        fi

        VCFFILE=${GUID}.stampy.vcf.gz
        if [ -f $VCFFILE ]; then
		echo 'vcf.gz output present'
	else
	        echo 'Mapping'
                # single-end mapping with Stampy
        	../../${STAMPY}stampy.py -g ../../mapto/tb_stampy -h ../../mapto/tb_stampy -M ${GUID}.fq.gz > stampy_output.sam 
                
	        ../../${SAMTOOLS}samtools view -t ../../../testdata/samtools_tab_file.txt  -o stampy_output.bam -b -S stampy_output.sam
	        rm -f stampy_output.sam
	        ../../${SAMTOOLS}samtools sort stampy_output.bam -t sort_tmpfile -o sorted_stampy_output.bam -m 8G
        	rm -f stampy_output.bam
        	../../${SAMTOOLS}samtools index sorted_stampy_output.bam

        	# can compute depth from INFO/AD tag
        	echo 'Variant calling'
	        ../../${SAMTOOLS}samtools mpileup -A -B -f ${REFFASTA} -Q25 -q30 -o40 -e20 -h100 -m2 -F0.002  -t INFO/AD -s -v  sorted_stampy_output.bam -o $VCFFILE
        	rm -f sorted_stampy_output.bam
       	 	rm -f sorted_stampy_output.bam.bai

	fi
    else
		echo Skipping $GUID as not a directory
    fi
    cd ..
done

exit 0



