# Map to reference genome

We provide an example script mapping fifty fastq files, obtained from a clinical *M. tuberculosis* sequencing operation, to the H37Rv reference genome.  

*1 Obtain the test data and software*    
The software required and test data can be obtained as [described](Prerequisites.md)

*2 Ensure mapper libraries are present    
Mapper libraries need to be generated in the mapto directory.  
The following should create both stampy and bowtie libraries.  
  
```
cd software
chmod +x prepare_to_map.sh
./prepare_to_map.sh
```

*3 Run mapping*       
A shell script is provided to do this.  

```
cd pipeline/testdata
chmod +x map_with_bowtie.sh

# recommend running with nohup, as it will take 30-60 minutes depending on the architecture
# note: the shell script check for output files; it will not rerun analyses
nohup ./map_with_bowtie.sh > map.out 2> map.err &
watch tail map.out

```




