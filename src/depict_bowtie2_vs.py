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
	rebuild_databases_if_present  = False,
	genbank_file_name = os.path.join("..", '..', "testdata", "NC_000962.3.gb")
	)


print("Depicting model")
am.depict_model()
print("Finished")

