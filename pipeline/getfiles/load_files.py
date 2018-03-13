# loads files from source data

import pandas as pd
import random
import os
    
inputfile = 'selected_meta.xlsx'
if os.path.exists(inputfile):
    selected_subset = pd.read_excel(inputfile)
    
else:
    # generate selected subset
    print('regenerating')
    inputfile = 'meta.xlsx'
    df = pd.read_excel(inputfile)
    selected_rows  = random.sample(df.index.tolist(),50)
    selected_subset =  df.loc[selected_rows]
    inputfile = 'selected_meta.xlsx'
    df.to_excel(inputfile)
 
    templates = {
        'vcf':
        'fastq':
        'bam':
        'kraken':  
    }


for i in df.index:
	sampleIds = list()
	missing = 0
	found =0
	failed =0
	for guid in guids:
		
		filename = "/home/dwyllie/dev/KRAKENREADER/kraken_data/{0}.kraken_report".format(guid)
		if not os.path.exists(filename):
			missing += 1
			cmd = "scp -P 8081 dwyllie@ndm.local@analysis1.mmmoxford.uk:/mnt/microbio/ndm-hicf/ogre/raw_input/{0}/QC/kraken.report /home/dwyllie/dev/KRAKENREADER/kraken_data".format(guid)
			os.system(cmd)
			try:
				os.rename("/home/dwyllie/dev/KRAKENREADER/kraken_data/kraken.report", "/home/dwyllie/dev/KRAKENREADER/kraken_data/{0}.kraken_report".format(guid))
				found +=1
			except OSError:		# doesn't exist
				print("Download failed. Path was /mnt/microbio/ndm-hicf/ogre/raw_input/{0}/QC/kraken.report".format(guid))
				f_missing.write("/mnt/microbio/ndm-hicf/ogre/raw_input/{0}/QC/kraken.report\n".format(guid))
				failed+=1
		else:
			found += 1
			#print("Found {0}".format(guid))
			#if found > 20:
			#	break

#print(selected_subset)
