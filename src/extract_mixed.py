#!/usr/bin/env python
""" extracts variation from regions within vcf files.

Is a command line wrapper for regionScan_from_genbank

Example usage:
python extract_mixed.py --help
python extract_mixed.py ../testdata/NC_000962.3.gb ../testdata/*v3.vcf.gz BaseCounts4 ../output

"""

# necessary libraries
import os
import glob
import argparse
from vcfScan import regionScan_from_genbank

if __name__ == '__main__':

    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('genbank_file_name', help='the path to a genbank format file containing details of the reference genome mapped to')
    parser.add_argument('input_file_pattern', help='the path to the input file(s).  Wildcards are accepted')
    parser.add_argument('infotag', help='the tag in the vcf INFO which contains the high-quality bases')
    parser.add_argument('outputdir', help='the directory for output')
    
    args = parser.parse_args()
    
    # check output directory exists
    if not os.path.exists(args.outputdir):
        raise ValueError("Output directory {0} does not exist".format(args.outputdir))
    
    # identify the genbank file of interest
    print("Extracting features from genbank file {0}".format(args.genbank_file_name))
    rs = regionScan_from_genbank(args.genbank_file_name, method = 'CDS', infotag=args.infotag)
    print("Feature extraction complete")	
    
    # export the extracted CDs to excel;
    output_excel = os.path.join(args.outputdir, 'regions.xlsx')
    rs.regions.to_excel(output_excel)
    print("Exported extracted features to {0}".format(output_excel)	)
    
    # identify files to process
    globpattern= os.path.join('..','testdata',"*v3.vcf.gz")
    inputfiles = glob.glob(args.input_file_pattern)
    print("Found {0} input files".format(len(inputfiles)))
    
    for inputfile in inputfiles:
        guid=os.path.basename(inputfile)[0:36]
        
        # test whether the file has already been parsed
        targetfile = os.path.join(args.outputdir,'{0}.txt'.format(guid))
        if os.path.exists(targetfile):
            print('{0} Target file exists {1}; skipped processing'.format(guid, targetfile))
    
        else:
            print("Examining file {0} ".format(guid))
            res = rs.parse(vcffile = inputfile, guid= guid)
            rs.region_stats.to_csv(targetfile)
            
