#!/bin/sh
# start in testdata;
cd ../testdata

# construct paths to executables relative to the pipeline directory
KRAKENDB=../../software/minikraken_20171101_4GB_dustmasked/
KRAKENSRCDIR=../../software/kraken-1.0/scripts/
REFFASTA=../../../testdata/NC_000962.3.fasta

echo 'Starting Kraken'

for dir in *
do
    GUID=${dir%*/}
    echo Analysing: ${GUID}
    if [ -d $GUID ]; then
        cd ${GUID}
        INPUTFILE=${GUID}.fq
	if [ -f ${INPUTFILE}.gz ]; then
		echo 'fastq.gz read file found; unzipping'
                gunzip ${INPUTFILE}.gz --to-stdout > ${INPUTFILE}
                ${KRAKENSRCDIR}kraken --threads 1 \
--preload --db $KRAKENDB \
--fastq-input $INPUTFILE |
kraken-filter --db \
$KRAKENDB --threshold 0.05 > ${GUID}.kraken.txt

kraken-report --db $KRAKENDB \
${GUID}.kraken.txt > ${GUID}.kraken.report.txt
                rm ${INPUTFILE}

	else
		echo $TARGETFILE not found
                exit 1
        fi
    else
		echo Skipping $GUID as not a directory
    fi
    cd ..
done

exit 0



