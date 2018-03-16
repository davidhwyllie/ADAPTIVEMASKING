#!/usr/bin/env python3
 
# necessary libraries
import os
from AdaptiveMasking import AdaptiveMasking

# assumes starting in pipeline/testdata
# fit model
print("Starting up")
am = AdaptiveMasking(
	analysis_name = 'bowtie',
	persistdir = os.path.join('..','output','bowtie'),
	genus_of_interest= 'Mycobacterium',
	rebuild_databases_if_present  = True,
	genbank_file_name = os.path.join("..", '..', "testdata", "NC_0103971.1.gb")
	)

print("Extracting region statistics from VCF files")
am.extract_model_input(vcf_inputpath=os.path.join('..','*.bowtie.vcf.gz')
					   
print("Reading in Kraken and vcf result files")
am.fit_model()

print("Depicting model")
am.read_model_input(kraken_inputpath = os.path.join('*','*.kraken.txt'),
				    region_inputpath = os.path.join('..','output','bowtie','*.regionstats.txt')
	)
