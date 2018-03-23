#!/usr/bin/env python3
 
# necessary libraries
import os
from AdaptiveMasking import AdaptiveMasking

# assumes starting in pipeline/testdata
# fit model
print("Starting up")
am = AdaptiveMasking(
	analysis_name = 'bowtie2_vs',
	persistdir = os.path.join('..','output','bowtie2_vs'),
	genus_of_interest= 'Mycobacterium',
	rebuild_databases_if_present  = True,
	genbank_file_name = os.path.join("..", '..', "testdata", "NC_000962.3.gb")
	)

print("Extracting region statistics from VCF files")
am.extract_model_input(vcf_inputpath=os.path.join('*','*.bowtie2_very_sensitive.vcf.gz'))
					   
print("Reading in Kraken and vcf result files")
am.read_model_input(kraken_inputpath = os.path.join('*','*.kraken.report.txt'),
				    region_inputpath = os.path.join('..','output','bowtie2_vs','*.regionstats.txt')
	)

print("Fitting model")
am.fit_model()
