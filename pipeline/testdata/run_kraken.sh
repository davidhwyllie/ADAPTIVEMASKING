#!/bin/sh
# run kraken across fastq.gz files
# start in testdata;
cd ../testdata

# construct paths to executables relative to the pipeline directory
KRAKENDB=../../software/minikraken_20171101_4GB_dustmasked/
KRAKENSRCDIR=../../software/kraken-1.0/scripts/
REFFASTA=../../../testdata/NC_000962.3.fasta

echo 'Starting Kraken'

for dir in $(ls */ -d)
do
    GUID=$(echo ${dir} | sed 's/.$//')
    echo Analysing: ${GUID}

    cd ${GUID}
    INPUTFILE=${GUID}.fq
    TARGETFILE=${GUID}.kraken.report.txt
    if [ -f ${TARGETFILE} ]; then
        echo 'Skipping; analysis already done'
    else
	if [ -f ${INPUTFILE}.gz ]; then
		echo 'fastq.gz read file found; unzipping'
               	gunzip ${INPUTFILE}.gz --to-stdout > ${INPUTFILE}
	        ${KRAKENSRCDIR}kraken --threads 1 \
--preload --db $KRAKENDB \
--fastq-input $INPUTFILE |
kraken-filter --db \
$KRAKENDB --threshold 0.05 > ${GUID}.kraken.txt

		kraken-report --db $KRAKENDB \
${GUID}.kraken.txt > ${TARGETFILE}
		# remove the temporary decompressed file
               	rm -f ${INPUTFILE}

	else
		echo ${INPUTFILE}.gz not found
               	exit 1
	fi
    fi
    cd ..
done

exit 0
