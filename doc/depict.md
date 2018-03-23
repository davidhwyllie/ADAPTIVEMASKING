# Depicting output

This section assumes that an analysis has been computed, as [described](model_maf.md).  The output from analyses is persisted, and can be viewed using an AdaptiveMasking object.  You instantiate this just as you would to do an analysis, but you *must* set rebuild_databases_if_present to False.
If you don't do this, any previous analysis will be deleted.


```python
from AdaptiveMasking import AdaptiveMasking

# read the results of a stored analysis
# critical: set rebuild_databases_if_present  = False,
# or the existing data will be overwritten
am = AdaptiveMasking(
	analysis_name = 'bowtie',
	persistdir = os.path.join('..','output','bowtie'),
	rebuild_databases_if_present  = False
)

```

You can now call methods which generate depictions, tables, and so on.

```
am.depict_model()

```

Example output is included, generated using a 250 sample test set and [four different mapping conditions](../example_output/example_output.md).
