# Model minor variant frequencies

A python class, *AdaptiveMasking*, performs many of the functions of the Adaptive Masking process.  This includes the modelling of minor variant frequencies as a function of non-target bacterial DNA quantities using Poisson regression models.  One model is constructed per gene.
The *fit_model()* method does this.  *read_model_input()* has to have been called beforehand.


The method used categorises non-target (i.e. non-Mycobacterial, in the example shown) bacterial DNA into a series of strata.  In the paper, we used strata bounded by 0,1,5,20 and 100%.  We chose these because they divided the samples into four approximately equal sizes, and elected to model the impact of the non-target bacterial DNA as a series of categories.
* The categories used can be controlled using the *categorisation_bins* option.
* The 'non-target' DNA is defined as being bacterial DNA (as assigned by Kraken) assigned to a genus or species which doesn't include *genus_of_interest*.

```
from AdaptiveMasking import AdaptiveMasking

am = AdaptiveMasking(
	analysis_name = 'bowtie',
	persistdir = os.path.join('..','output','bowtie'),
	genus_of_interest= 'Mycobacterium',
	rebuild_databases_if_present  = True,
	genbank_file_name = os.path.join("..", '..', "testdata", "NC_0103971.1.gb"),
    categorisation_bins = [0, 0.01, 0.05, 0.2, 1]
	)

```

If not done already, region statistics can be extracted
```
print("Extracting region statistics from VCF files")
am.extract_model_input(vcf_inputpath=os.path.join('*','*.bowtie.vcf.gz'))
```

And Kraken and vcf region information read in
```
print("Reading in Kraken and vcf result files")
am.read_model_input(kraken_inputpath = os.path.join('*','*.kraken.report.txt'),
				    region_inputpath = os.path.join('..','output','bowtie','*.regionstats.txt')
	)
```

After this, the models can be fitted:
```
print("Fitting model")
am.fit_model()
```
Depiction is described [here](depict.md).

