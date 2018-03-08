#!/usr/bin/env python3
 
# necessary libraries
import os
import pandas as pd
import glob
import numpy as np
import logging
from KrakenReportReader import KrakenReportReader
import statsmodels.api as sm

# define the input variables
class AdaptiveMasking():
	""" computes minor variant frequencies across genomic regions,
	and their relationship to non-target DNA present in the sample """
	
	def __init__(self,
				 analysis_name,
				 persistdir,
				 rebuild_databases_if_present = True,
				 genus_of_interest = 'Mycobacterium',
				 categorisation_bins = [0, 0.01, 0.05, 0.2, 1],
				 model_min_depth = 0
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
			Returns:
				None
		"""
		
		self.analysis_name= analysis_name	
		self.persistdir = persistdir
		self.genus_of_interest = genus_of_interest
		self.categorisation_bins = categorisation_bins
		self.model_min_depth = model_min_depth
			
		# now recover all the data from the genetic parsing.
		# store in an hdf5 data structure.
		self.hdf_file = os.path.join(self.persistdir,'{0}.h5'.format(analysis_name))
		if os.path.exists(self.hdf_file) and rebuild_databases_if_present:
			os.unlink(self.hdf_file)
			
		return None	
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
			sampleId = os.path.basename(inputfile)[0:36]
			result = krr.simplify(inputfile = inputfile, guid = sampleId)
			results.append(result)
			
			nRead+=1
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
			guid = df['sampleId'].unique()[0]
		
			if guid in kraken_sampleIds:		# if there is a kraken report
				# read into database
		
				df.to_hdf(self.hdf_file, 'model_input',
						  mode='a',
						  format='table',
						  append= True,
						  data_columns= ['sampleId','gene'])
				nRead+=1
				if nRead % 100 == 0:
					logging.info("Added {0}".format(nRead))
					
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
			
		# recovery kraken_pseudo, which contains the explanatory data
		kraken_pseudo = pd.read_hdf(self.hdf_file, 'explan')
		
		# recover all genes in the database
		res = pd.read_hdf(self.hdf_file, 'model_input', columns=['gene'])
		roi_names = res['gene'].unique()
		n_rois = 0
		for roi_name in roi_names:
			
			## construct a statistical model for each roi_name

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
				print("Model fitting failed")
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
			if n_rois > 50:
				break
		
		estimated_mixtures.to_hdf(self.hdf_file, 'model_output',
						  mode='a',
						  format='table',
						  append= False)
		self.coefficients = estimated_mixtures
		return estimated_mixtures
		
am = AdaptiveMasking(
	analysis_name = 'test1',
	persistdir = os.path.join('..','modelling','tmp'),
	genus_of_interest= 'Mycobacterium'
					)
am.read_model_input(kraken_inputpath = os.path.join('..','modelling','kraken_reports','*.kraken_report'),
				    region_inputpath = os.path.join('..','modelling','region_reports','*.tsv'),
					max_files_analysed = 50
					)
coefficients = am.fit_model()
print(coefficients)
exit()



print(estimated_mixtures)
