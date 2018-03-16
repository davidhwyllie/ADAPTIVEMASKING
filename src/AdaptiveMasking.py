#!/usr/bin/env python3
 
# necessary libraries
import os
import pandas as pd
import glob
import numpy as np
import logging

import statsmodels.api as sm
from bokeh.layouts import gridplot
from bokeh.plotting import figure, show, output_file
from bokeh.io import output_file, show
from bokeh.layouts import row, column, gridplot, widgetbox
from bokeh.models import  NumberFormatter
from bokeh.models import  ColumnDataSource,\
                          Circle, \
                          HoverTool, SaveTool, ResetTool, BoxZoomTool,\
                          ZoomInTool, ZoomOutTool,  PanTool, TapTool,\
                          WheelZoomTool, BoxSelectTool, PolySelectTool,\
                          Span, Arrow, OpenHead
from bokeh.models.widgets import DataTable, DateFormatter, \
                          TableColumn, CheckboxGroup, RadioButtonGroup, \
                          TextInput, Panel, Tabs
from KrakenReportReader import KrakenReportReader
from vcfScan import regionScan_from_genbank

class AdaptiveMasking():
	""" computes minor variant frequencies across genomic regions,
	and their relationship to non-target DNA present in the sample """
	
	def __init__(self,
				 analysis_name,
				 persistdir,
				 rebuild_databases_if_present = True,
				 genus_of_interest = 'Mycobacterium',
				 categorisation_bins = [0, 0.01, 0.05, 0.2, 1],
				 model_min_depth = 0,
				 genbank_file_name = None
				 ):
		""" prepare to compute minor variant frequencies across genomic regions,
			and their relationship to non-target DNA present in the sample.
			
			Arguments:
				analysis_name: an alphanumeric designation for the analysis.
								Used to name the hdf file in which the analysis is persisted.
				rebuild_databases_if_present: (default True)
								if False, will not recreate the data used for depiction.
				persistdir: a writable directory used for storing an hdf5 database containing the analysis.
				genus_of_interest:  a string, such as 'Mycobacterium' or 'Mycobacterium tuberculosis' the presence of which in the
							scientific name of an organism identifies it as being the 'expected result', as opposed to being 'non-expected' bacterial DNA
				categorisation_bins:  a list containing the limits of the categories of non-expected bacterial DNA which are modelled
				model_min_depth: do not fit models if observations have depth less than some cutoff.  minimum zero.
				genbank_file_name: the genbank file to parse to obtain the regions of interest.  Must be specified if rebuild_database_if_present = True, or the database doesn't exist.
			Returns:
				None
		"""
		
		# configure logging;
		logger = logging.getLogger()
		logger.level = logging.INFO
		
		# persist parameters as object properties;
		self.analysis_name= analysis_name	
		self.persistdir = persistdir
		self.genus_of_interest = genus_of_interest
		self.categorisation_bins = categorisation_bins
		self.model_min_depth = model_min_depth
			
		# configure an hdf5 data structure.
		self.hdf_file = os.path.join(self.persistdir,'{0}.h5'.format(analysis_name))
		if os.path.exists(self.hdf_file) and rebuild_databases_if_present:
			os.unlink(self.hdf_file)
		if not os.path.exists(self.hdf_file):		# we're making in a new one
			
			# create a new regionScan_from_genbank object.
			# run checks on genbank file name
			if genbank_file_name is None:
				raise ValueError("genbank_file_name not specified, but is required for new object creation")
			if not os.path.exists(genbank_file_name):
				raise FileExistsError("genbank_file_name specifies non-existent file {0}".format(genbank_file_name))

			# generate a regionScan_from_genbank object for the first time.
			self.rs = regionScan_from_genbank(genbank_file_name, method = 'CDS')

			# persist self.rs.regions, which are the regions to be examined.
			self.rs.regions.to_hdf(self.hdf_file, 'regions',
							  mode='w',
							  format='table',
							  append= False,
							  data_columns= ['roi_name','mid','start_pos','end_pos'])

		else:   # we generate a regionScan_from_genbank object based on stored data
			
			# recover self.regions from hdf file.
			# TODO: make a new HDF
			regions =  pd.read_hdf(self.hdf_file, 'regions')
			self.rs = regionScan_from_genbank(genbank_file_name=regions)

		return None	
	def extract_model_input(self,
						 vcf_inputpath,
						 max_files_analysed = 1e8
						):
		""" analyses vcf files found on vcf_inputpath and writes the region summaries to persistdir.
		
		wraps regionScan_from_genbank.
		
		Arguments:
				region_inputpath: a path which when globbed yields vcf files for analysis.
				max_files_analysed: the maximum number of files which will be analysed.  Useful only for debugging/demos.
				
		Returns:
				none
				Output is written to persistdir.
		"""
		
		# identify files to process
		inputfiles = glob.glob(vcf_inputpath)
		logging.info("Found {0} input files".format(len(inputfiles)))
		nRead=0
		
		for inputfile in inputfiles:
			guid=os.path.basename(inputfile)[0:36]
			
			# test whether the file has already been parsed
			targetfile = os.path.join(self.persistdir,'{0}.regionstats.txt'.format(guid))
			if os.path.exists(targetfile):
				print('{0} Target file exists {1}; skipped processing'.format(guid, targetfile))
		
			else:
				logging.info("Examining file {0} ".format(guid))
				res = self.rs.parse(vcffile = inputfile, guid= guid)
				self.rs.region_stats.to_csv(targetfile)

				nRead+=1
				if nRead % 5== 0:
					logging.info("Read {0} per-gene report files".format(nRead))
					
				if nRead > max_files_analysed:
					logging.warning("Stopped parsing as max_files_analysed was reached")
					break
	def read_model_input(self,
						 kraken_inputpath,
						 region_inputpath = os.path.join('..','modelling','region_reports','*.tsv'),
						 max_files_analysed = 1e8
						):
		""" reads in the inputs to the statistical model.

		Arguments:
				kraken_inputpath: a path which when globbed yields kraken reports
				region_inputpath: a path which when globbed yields region reports
		Returns:
				none
				Output is written to the hdf5 database at self.hdf_file
		"""

		# make KrakenReader object, identifying the genus Mycobacterium.
		krr = KrakenReportReader(genus_of_interest = self.genus_of_interest)

		results = []
		nRead = 0
		logging.info("Summarising Kraken data for analysis")
		for inputfile in glob.glob(kraken_inputpath):
			print(inputfile)
			sampleId = os.path.basename(inputfile)[0:36]
			result = krr.simplify(inputfile = inputfile, guid = sampleId)
			results.append(result)

			nRead+=1
			if nRead % 100 == 0:
				logging.info("Read {0} Kraken files".format(nRead))
			if nRead > max_files_analysed:
				logging.warning("Stopped reading Kraken data as max_files_analysed was reached")
				break

		# generate a data frame comprising the Kraken summary
		kraken = pd.DataFrame.from_dict(results)

		# if there are no reads, then we cannot analyse such samples.  we drop them.
		kraken = kraken[kraken['nReads']>0]
		logging.info("Excluded results with zero reads.  {0} remain.".format(len(kraken.index)))


		# categorise the proportion of bacterial reads which are not of interest.
		# and create pseudovariables while dropping the first column which is the reference category.
		for i in kraken.index:
			kraken.loc[i, 'prop_bug'] = kraken.loc[i,'B']/kraken.loc[i,'nReads']

		kraken['prop_bug_cat'] =  pd.cut(kraken['prop_bug'], bins=self.categorisation_bins, include_lowest=True)
		kraken_pseudo = pd.crosstab(kraken['sampleId'], kraken['prop_bug_cat'])
		kraken_pseudo.columns = kraken_pseudo.columns.astype(str)

		# drop the first column: it is reference
		kraken_pseudo = kraken_pseudo.drop(axis=1, columns = kraken_pseudo.columns.tolist()[0])
		kraken_pseudo.to_hdf(self.hdf_file, 'explan',
						  mode='a',
						  format='table',
						  append= True,
						  data_columns= ['sampleId','gene'])
		logging.info("Recovered {0} kraken results".format(len(kraken.index)))

		# now recover all the data from the genetic parsing.
		# store in an hdf5 data structure.
		results = []
		nRead = 0
		kraken_sampleIds = set(kraken['sampleId'])

		inputfiles = glob.glob(region_inputpath)
		for inputfile in inputfiles:
			df = pd.read_csv(inputfile, sep='\t')
			guid = basename(inputfile)[0:38]
			df['sampleId']= guid
		
			if guid in kraken_sampleIds:		# if there is a kraken report
				# read into database
				try:
					df.to_hdf(self.hdf_file, 'model_input',
							  mode='a',
							  format='table',
							  append= True,
							  data_columns= ['sampleId','gene'])
				except ValueError:
					# occurs if the data file is malformed or missing
					logging.warn("Data file could not be loaded; malformed {0}".format(inputfile))
					print(df)
					
				nRead+=1
				if nRead % 100 == 0:
					logging.info("Read {0} per-gene report files".format(nRead))
					
				if nRead > max_files_analysed:
					logging.warning("Stopped reading Kraken data as max_files_analysed was reached")
					break

	def fit_model(self):
		""" fits poisson regression model relating the
		minor variant frequency to the amount of non-target bacterial DNA
		
		Arguments:
			None
			
		Returns:
			a data frame containing coefficients from the modelling.
			This is also stored as self.coefficients, and can be exported to excel, csv etc.
			"""
		
		logging.info("Fitting models; reading Kraken data")	
		# recovery kraken_pseudo, which contains the explanatory data
		try:
			kraken_pseudo = pd.read_hdf(self.hdf_file, 'explan')
		except KeyError:
			raise KeyError("Tried and failed to read Kraken data.  You must load Kraken data with .read_model_input() before fitting models.")
		
		# recover all genes in the database
		logging.info("Fitting models; reading gene names")	
		res = pd.read_hdf(self.hdf_file, 'model_input', columns=['gene'])
		roi_names = res['gene'].unique()
		n_rois = 0
		for roi_name in roi_names:		# ['B55','rrs','rrl']:
			
			## construct a statistical model for each roi_name
			logging.info("Modelling {0}".format(roi_name))
			
			# construct data frame for modelling
			model_input = pd.read_hdf(self.hdf_file, 'model_input',  where='gene=="{0}"'.format(roi_name))
			model_input = model_input[['sampleId','total_depth','total_nonmajor_depth']]
			model_input = model_input.set_index('sampleId',drop=True)
			model_input = kraken_pseudo.merge(model_input, how='inner', left_index=True, right_index=True)	
			explan_cols = kraken_pseudo.columns.tolist()
			model_input = model_input[model_input['total_depth']> self.model_min_depth]
			model_input['offset'] = [np.log(int(x)) for x in model_input['total_depth']]
			explan_vars = model_input[explan_cols]
			explan_vars = sm.add_constant(explan_vars, prepend=True)
			
			if sum(model_input['total_nonmajor_depth']) ==0:
				# then initial setting of mu fails.  looks to me like regression of this bug:
				# https://bugs.launchpad.net/statsmodels/+bug/603306
				# pending fix in statsmodels, in this situation we set the total_nonmajor_depth to 1 for one item.
				model_input.loc[model_input.index[0], 'total_nonmajor_depth'] = 1
		
			## fit model
			# total_nonmajor_depth ~ prop_bug_catcat + offset(log(total_depth)
		
			try:
				poisson_model = sm.GLM(model_input['total_nonmajor_depth'],
									   explan_vars,
									   offset = model_input['offset'],
									   family = sm.families.Poisson()
									)
				poisson_results  = poisson_model.fit()
				model_fit_succeeded = True
				
			except sm.tools.sm_exceptions.PerfectSeparationError:
				# there isn't enough data and colinearity prevents Hessian computation;
				model_fit_succeeded = False
			
			if model_fit_succeeded == True:	
				res = pd.DataFrame({
					'Fitting_succeeded':True,
					'Estimate':poisson_results.params,
					'p_value':poisson_results.pvalues,
					'lower_ci':poisson_results.conf_int()[0],
					'upper_ci':poisson_results.conf_int()[1]}
					)
				res['nObs'] = poisson_results.nobs
				res['roi_name'] = roi_name
				res['parameter']=res.index
				
			else:
				logging.WARN("Model fitting failed, skipping region {0}".format(roi_name))
				model_input.to_hdf(hdf_file, 'modelling_failed_input',
						  mode='a',
						  format='table',
						  append= True)
				res = pd.DataFrame({
					'Fitting_succeeded':False,
					'Estimate':None,
					'p_value':None,
					'lower_ci':None,
					'upper_ci':None,
					'nObs':len(model_input.index),
					'roi_name':roi_name,
					'parameter':'const'},
					index = [roi_name]
					)
				
			n_rois +=1
			
			if n_rois == 1:
				estimated_mixtures = res
			else:
				estimated_mixtures = estimated_mixtures.append(res, ignore_index=True)
				print("Fitting {0} #{1}".format(roi_name, n_rois))
			
			# debug
			#if n_rois > 50:
			#	break
		
		estimated_mixtures.to_hdf(self.hdf_file, 'model_output',
						  mode='a',
						  format='table',
						  append= False)
		self.coefficients = estimated_mixtures
		return estimated_mixtures
	def depict_model(self,
					 exclude_estimates_over = 10000):
		
		""" depicts the coefficients fitted.
		
			Arguments:
				
				exclude_estimates_over:  will not plot parameter estimates more than this number.
				In this setting, the background rate of mixtures is about 0.001 and the
				'estimate' referred to is an incidence rate ratio relative to the background rate;
				Observed estimates for highly affected genes are in the range 5-15.
				Useful for exclusion of rare parameter estimates which are essentially infinite due to non-convergence.
			
			"""
			
		
		# read region annotation from genbank entry

		try:
			self.coefficients = pd.read_hdf(self.hdf_file, 'model_output')
		except KeyError:
			raise KeyError("Attempted to read model output from file {0} but no model output is present.  you need to fit the model using .fit_model() before attempting to depict the output".format(self.hdf_file))

		# get the parameters modelled
		param_values = self.coefficients['parameter'].unique()
		highest_test_cat = max(set(param_values)-set(['const']))

		# compute the IRR/ fold change (for the coefficients)
		self.coefficients['exp_Estimate'] = [np.exp(x) for x in self.coefficients['Estimate']]
		self.coefficients['log10_Estimate'] = [2.303*x for x in self.coefficients['Estimate']]

		df1 = self.coefficients.query("parameter=='{0}'".format(highest_test_cat))
		df2 = self.coefficients.query("parameter=='{0}'".format('const'))
		df1 = df1[['exp_Estimate','roi_name']]
		df2 = df2[['exp_Estimate','roi_name']]
		df1 = df1.rename({'exp_Estimate':'hi'}, axis=1)
		df2 = df2.rename({'exp_Estimate':'const'}, axis=1)
		bivar = df1.merge(df2, how = 'inner', on= 'roi_name')
		bivar = bivar.merge(self.rs.regions, how = 'inner', on= 'roi_name')

		bivar_source= ColumnDataSource(bivar)
		bivar_tools = "save,reset,box_zoom,zoom_in,zoom_out,pan,tap,box_select"
		
		bp = figure(plot_width=600, plot_height=400,tools=bivar_tools,toolbar_location='right')
		bp.title.text = "Mixtures in the presence vs. absence of unexpected bacterial DNA"
		r0 = bp.circle(source = bivar_source,  x='const', y='hi', size= 2, alpha=0.5)
		bp.xaxis.axis_label = 'Est. Mixture when reference category amount of unexpected bacterial DNA'
		bp.yaxis.axis_label = 'IRR mixture estimate with {0} extr. DNA'.format(highest_test_cat)
		
		g1 = figure(plot_width=600, plot_height=200,tools=bivar_tools,toolbar_location='right')
		g1.title.text = "Gene positions"
		r1 = g1.circle(source = bivar_source,  x='mid', y='hi', size= 2, alpha=0.5)
		g1.xaxis.axis_label = 'Genome position'
		g1.yaxis.axis_label = 'IRR mix. {0} extr. DNA'.format(highest_test_cat)

		g2 = figure(plot_width=600, plot_height=200,tools=bivar_tools,toolbar_location='right')
		g2.title.text = "Gene positions"
		r2 = g2.circle(source = bivar_source,  x='mid', y='const', size= 2, alpha=0.5)
		g2.xaxis.axis_label = 'Genome position'
		g2.yaxis.axis_label = 'Est. mix with ref cat. extr. DNA'.format(highest_test_cat)
		 
		hover1=HoverTool(renderers=[r0, r1,r2])
		hover1.tooltips= [("roi_name", "@roi_name"),("% mix ref. cat", "@const"),("IRR mix {0}".format(highest_test_cat),"@hi")]
		bp.add_tools(hover1)
		g1.add_tools(hover1)
		g2.add_tools(hover1)				

		# define columnar data source
		          # define columnar data source
		columns = [
				  TableColumn(field="roi_name", title="Region analysed"),
				  TableColumn(field="const", title="Mix with ref. cat extra DNA", formatter = NumberFormatter(format ='0.00000')),
				  TableColumn(field='hi', title='IRR mix with {0} extra DNA'.format(highest_test_cat), formatter = NumberFormatter(format = '0.00')),
				  TableColumn(field='start_pos', title='region start', formatter = NumberFormatter(format = '0')),
				  TableColumn(field='length', title='region length', formatter = NumberFormatter(format = '0'))	
				  ]
		data_table = DataTable(source=bivar_source,
							   columns=columns,
							   width=600,
							   height=200,
							   scroll_to_selection= True,
							   sortable = True,
							   reorderable = True,
							   editable = False,
							   fit_columns = True)
		
		lc = column ( [row([bp]), row([g1]), row([g2])	] )	# row([data_table]), 
		content = Panel(child = lc, title='Bivariate')
			
		# add the plots to a tab
		tab_content = []
		tab_content.append(content)
		tls = Tabs (tabs = tab_content)

		# compute histograms for each parameter
		for param_value in param_values:
	
			if param_value == 'const':
				x_axis_label = 'Minor variant frequency with reference category unexpected bacterial DNA'
				tab_name = 'Low extr. DNA'
			else:
				x_axis_label = 'IRR if unexpected bacterial DNA in range {0} relative to ref. cat'.format(param_value)
				tab_name = 'Extr. DNA '+param_value
			these_coeffs = self.coefficients.query("parameter == '{0}'".format(param_value))
			
			
			these_coeffs = these_coeffs[these_coeffs['Fitting_succeeded']==True]

			# sanity check: these_coeffs must be < exclude_estimates_over;
			# if not, estimates are likely to be unstable
			# only one coefficient is excluded in the test set with a cutoff of 50
			these_coeffs = these_coeffs[these_coeffs['exp_Estimate']< exclude_estimates_over]

			# construct histogram and cusum
			hist,edges = np.histogram(these_coeffs['exp_Estimate'], bins=300)
			cumdist = [0]*len(hist)
			for i in range(len(cumdist)-1):
				cumdist[i+1]= hist[i]+cumdist[i]
			for i in range(len(cumdist)):
				cumdist[i]= cumdist[i]/sum(hist)
	
			# depict
			hist_tools = "save,reset,box_zoom,zoom_in,zoom_out,pan"
			p1 = figure(title="Estimated Minor variant frequency across regions",tools=hist_tools,
						background_fill_color="#E8DDCB", plot_height=300)
			p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
			fill_color="#036564", line_color="#033649")
			p1.legend.location = "center_right"
			p1.legend.background_fill_color = "darkgrey"
			p1.xaxis.axis_label = x_axis_label
			p1.yaxis.axis_label = 'Number of regions'
			
			if not param_value=='const':
				vline = Span(location=1, dimension='height', line_color='red', line_width=1)
				p1.renderers.extend([vline])
			
			p2 = figure(title="Cumulative Estimated Minor variant frequency across regions",tools=hist_tools,
						background_fill_color="#E8DDCB", plot_height=300)
			p2.quad(top=cumdist, bottom=0, left=edges[:-1], right=edges[1:],
			fill_color="#036564", line_color="#033649")
			p2.legend.location = "center_right"
			p2.legend.background_fill_color = "darkgrey"
			p2.xaxis.axis_label = x_axis_label
			p2.yaxis.axis_label = 'Cumulative proportion of regions'
	
			lc = column ( [row([p1]), row([p2])])
			content = Panel(child = lc, title=tab_name)
			
			# add the plots to a tab
			tab_content.append(content)
		tls = Tabs (tabs = tab_content)

		layout = tls
		show(layout) 

