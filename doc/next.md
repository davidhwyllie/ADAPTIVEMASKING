# Next

Having identified regions which you believe need masking due to possible low-confidence calling in these region(s), you can take the following actions:

1) Prior to phylogenetic tree construction
If you have a sufficiently small number of samples that you can generate phylogenetic trees, 'mask' any consensus bases called in the regions by replacing them with N.
This tells the phylogenetic tree generation software that no information is contained within these regions.

2) Find similar sequences using SNV based sequence comparisons, while ignoring low regions which you think (based on analysis of modelled minor variant frequencies) should be masked
We have build a [tool](https://github.com/davidhwyllie/findNeighbour2) to do this.  This operates as a server, and exposes a RESTful API for sequence addition and data querying.
When the server is configured, an [EXCLUDEFILE](https://github.com/davidhwyllie/findNeighbour2/blob/master/doc/HowToTest.md) containing bases to be excluded can be provided; excluded bases are ignored, whether or not 
they are set to 'N'.


