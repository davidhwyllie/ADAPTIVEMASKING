# selects input file for demonstration pipeline from the MMM computational environment.

import pandas as pd
import random
import os
import shutil

inputfile = 'meta.xlsx'
df = pd.read_excel(inputfile)

templates = {
    'vcf.gz':'/mnt/microbio/ndm-hicf/ogre/pipeline_output/{0}/MAPPING/2e6b7bc7-f52c-4649-8538-c984ab3894bb_R00000039/STD/basecalls/{0}_v3.vcf.gz',
    'bam':'/mnt/microbio/ndm-hicf/ogre/raw_input/{0}/in.bam',  
    'kraken':'/mnt/microbio/ndm-hicf/ogre/raw_input/{0}/QC/kraken.report'  
}

guids_found = 0
selected_rows = []
for i in df.index:
    
    guid = df.loc[i, 'sequenceGuid']
    nFound = 0
    export = []
    destdir = os.path.join("..", 'testdata',guid)
    if  not os.path.exists(destdir):
        for filetype in templates.keys():
            expected_file = templates[filetype].format(guid)
            target_file = os.path.join(destdir, '{0}.{1}'.format(guid,filetype))
            export.append((expected_file, target_file))
            #print(guid, expected_file, os.path.exists(expected_file))
            if os.path.exists(expected_file):
                nFound +=1
                
        if nFound == 3:
            print(guid)
            os.mkdir(destdir)
            print("Created", destdir)
            for (expected_file, target_file) in export:
                shutil.copy(expected_file, target_file)
            guids_found +=1
            selected_rows.append(i)
    if guids_found == 50:
        break
    

selected_subset =  df.loc[selected_rows]
inputfile = 'selected_meta.xlsx'
df.to_excel(inputfile)
