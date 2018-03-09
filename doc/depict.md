# Depicting output

This section assumes that an analysis has been computed, as [described](model_maf.md).  The output from analyses is persisted, and can be viewed using an AdaptiveMasking object.  You instantiate this just as you would to do an analysis, but you *must* set rebuild_databases_if_present to False.
If you don't do this, any previous analysis will be deleted.


```python
from AdaptiveMasking import AdaptiveMasking

# read the results of a stored analysis
# critical: set rebuild_databases_if_present  = False,
# or the existing data will be overwritten
am = AdaptiveMasking(
	analysis_name = 'test1',
	persistdir = os.path.join('..','modelling','tmp'),
	genus_of_interest= 'Mycobacterium',
	rebuild_databases_if_present  = False
					)

```

You can now call methods which generate depictions, tables, and so on.

```
am.depict_model()

```
