# Extract minor variant frequencies from VCF files

This is done using the regionScan_from_genbank class.
The code examines a VCF file, of which an example is below.
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EXTRA
R00000039	1	.	T	G	38.10	K0.90;z	ABQ4=39.250;BCALL=N;BaseCounts=0,0,7,81;BaseCounts4=0,0,5,27;DM4=-3.103;DM4L=-7.856;DP=88;DP4=24,3,5,0;DPT4L=-81.656;DZ4=-2.480;DZ4L=-6.927;GC=60.780;MQ=59;MQ4=60;PCALL4=0.000;PCONS4=1.000;SBR=0	GT:DP	0/1:32
R00000039	2	.	T	G	32.10	K0.90;z	ABQ4=38.824;BCALL=N;BaseCounts=0,0,8,81;BaseCounts4=0,0,5,29;DM4=-3.013;DM4L=-7.947;DP=89;DP4=26,3,5,0;DPT4L=-82.759;DZ4=-2.408;DZ4L=-7.005;GC=59.620;MQ=59;MQ4=60;PCALL4=0.000;PCONS4=1.000;SBR=0	GT:DP	0/1:34
R00000039	3	.	G	T	14.20	K0.90;S25;z	ABQ4=38.750;BCALL=N;BaseCounts=0,1,83,5;BaseCounts4=0,0,28,4;DM4=-3.103;DM4L=-7.944;DP=90;DP4=25,3,4,0;DPT4L=-82.726;DZ4=-2.480;DZ4L=-7.003;GC=60.380;MQ=59;MQ4=60;PCALL4=0.000;PCONS4=1.000;SBR=0	GT:DP	0/1:32
```

It expects a field in the INFO section which contains the numbers of A,C,G,T to process.  
In the above example, the BaseCounts and BaseCounts4 tags are examples.

A working example, using the project's test data, is below.  The code is supplied as working_example_of_maf_extraction.py

```python
#!/usr/bin/env python
""" example of use of regionScan """

# necessary libraries
import os
import glob
from vcfScan import regionScan_from_genbank


# identify the genbank file of interest
genbank_file_name = os.path.join("..", "testdata", "NC_000962.3.gb")
print("Extracting features from genbank file {0}".format(genbank_file_name))
rs = regionScan_from_genbank(genbank_file_name, method = 'CDS', infotag='BaseCounts4')
print("Feature extraction complete")	

```

This stores the regions of interest in the regionScan_from_genbank object.
The regions generated can be exported:

```python
# export the extracted rois to excel;
output_excel = os.path.join('..','output', 'regions.xlsx')

# note that rs.regions is a pandas dataframe, and pandas methods can be called on it;
rs.regions.to_excel(output_excel)
print("Exported extracted features to {0}".format(output_excel)	)

```

The regions extracted from the genbank file look something like this:

```
start_pos	end_pos	name
1	1524	Rv0001
1525	2051	near_Rv0002
2052	3260	Rv0002
3261	3279	near_Rv0003
3280	4437	Rv0003
4438	4433	near_Rv0004

```

Now, we are ready to parse the input vcf files.
An example of how to do this follows, outputting the summary data into a series of .csv files:

```python
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
        
        # Optional: If you want subsequent random access to subsections of the
        # regions examined, please see the .persist() method in the vcfScan class.
        # This will persist per-base information in the regions examined (which, in the
        # above case, includes the whole genome) to an indexed hdf5 format.

```


The output per region looks like this:
```
roi_name,mean_depth,min_depth,max_depth,start,stop,length,mean_maf,total_depth,total_nonmajor_depth
Rv0001,73.99540682414698,21,126,1,1524,1524,0.0005167313777354538,112769,60
Rv0002,41.789909015715466,10,91,2052,3260,1209,0.00108890037274931,50524,47
Rv0003,64.38082901554404,33,94,3280,4437,1158,0.0006818405450617134,74553,52
Rv0004,55.28900709219858,32,82,4434,4997,564,0.0009366414816737341,31183,28
Rv0005,70.75493096646943,41,96,5240,7267,2028,0.0006573041881978442,143491,96
Rv0006,63.11601112435439,41,96,7302,9818,2517,0.0008032353616680821,158863,129
Rv0007,42.90819672131148,17,75,9914,10828,915,0.0005097938211881642,39261,19
Rv0008c,44.550228310502284,26,63,11874,12311,438,0.0004895420366389829,19513,9
```

The output order is sorted by roi_name.
Field names are as follows:

* roi_name the name of the region of interest
* mean_depth the arithmetic mean depth, defined as the sum of A,C,G,T in the INFO infotag section
* min_depth  the arithmetic min depth, defined as the sum of A,C,G,T in the INFO infotag section
* max_depth  the arithmetic max depth, defined as the sum of A,C,G,T in the INFO infotag section
* start the 1-indexed start of the roi_name in the reference genome, inclusive
* stop  the 1-indexed end of the roi_name in the reference genome, inclusive
* length the length of the region examined in nucleotides
* mean_maf  maf is defined as the second most common nucleotide in the INFO infotag section/ the total depth.  mean_maf is the sum of maf / length.
* total_depth the sum of the per-base depths across the region of interest
* total_nonmajor_depth the sum of the per-base depths, except for the most common base, across the region of interest.

In the approach taken in the linked paper,  
total_nonmajor_depth is modelled as a function non-Mycobacterial bacterial DNA concentrations, using total_depth as an offset.


A freestanding script suitable for command line use, extract_mixed.py, is also available.
Please see the documentation in the file.

The following command would do the same as the example above:

```
python extract_mixed.py ../testdata/NC_000962.3.gb ../testdata/*v3.vcf.gz BaseCounts4 ../output

```
