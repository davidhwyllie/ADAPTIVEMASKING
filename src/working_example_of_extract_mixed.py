#!/usr/bin/env python
""" example of use of regionScan """

# necessary libraries
import os
import glob
from vcfScan import regionScan_from_genbank


# identify the genbank file of interest
genbank_file_name = os.path.join("..", "testdata", "NC_000962.3.gb")
print("Extracting features from genbank file {0}".format(genbank_file_name))
rs = regionScan_from_genbank(genbank_file_name, method = 'gene', infotag= 'BaseCounts4')
print("Feature extraction complete")	

# export the extracted CDs to excel;
output_excel = os.path.join('..','output', 'regions.xlsx')

# note that rs.regions is a pandas dataframe, and pandas methods can be called on it;
rs.regions.to_excel(output_excel)
print("Exported extracted features to {0}".format(output_excel)	)

# define where the output is to go.
outputdir = os.path.join('..','output')

# identify files to process
globpattern= os.path.join('..','testdata',"*v3.vcf.gz")
inputfiles = glob.glob(globpattern)
print("Found {0} input files".format(len(inputfiles)))


for inputfile in inputfiles:
    guid=os.path.basename(inputfile)[0:36]
    
    # test whether the file has already been parsed
    targetfile = os.path.join(outputdir,'{0}.txt'.format(guid))
    if os.path.exists(targetfile):
        
        print('exists {0}'.format(guid))

    else:
        print("Examining file {0} ".format(guid))
        res = rs.parse(vcffile = inputfile, guid= guid)
        rs.region_stats.to_csv(targetfile)
        
