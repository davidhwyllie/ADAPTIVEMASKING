#!/usr/bin/env python3
"""
classes to parse VCF files extracting high-quality bases from a pre-defined set of positions.

Applications include detection of mixtures of TB genomes at lineage defining sites, or sites varying within a genome.
These classes do not create a persistent, indexed extract of a subset of data and do not allow
rapid selection of regions of mixed base calls.  Please see the TBMIX class if you need to do this.

To modify to read BCF files, modify the ._parse library to use pysam module to for file access.

Note that in general, some methods within these classes require a sequence identifier to be provided.
These are referred to as 'guids' and we hav always used globally unique identifiers (guids) as sequence identifiers.
However, the requirement to actually use a guid is not absolute, and is not enforced.
Other character strings identifying a sequence should also work, subject to:
* not being more than 36 characters
* being valid as part of a file name
"""

import csv
import os
import sys
import gzip
import shutil
import glob
import json
import hashlib
from collections import deque
from scipy import stats
import numpy  as np
import time
import math
import unittest
import logging
import warnings
import pandas as pd
import copy
import datetime
import warnings
import tables
from Bio import SeqIO, SeqRecord, SeqFeature


class vcfScan():
	""" parses a VCF file, extracting high quality calls at pre-defined positions.
	
		Note that at present this doesn't support BCF files.
		To enable this, we'd need to modify the ._parse function
		to use pysam https://pypi.python.org/pypi/pysam.
		
		
	"""
	
	def __init__(self, 
	                expectedErrorRate = 0.001,
	                infotag = 'BaseCounts4'):
		""" creates a vcfScan object. 
		
		Args:
		    expectedErrorRate: a floating point number indicating the expected per-base error rate.
		    infotag: a string indicating the vcf/bcf INFO tag from which to extract the four high quality depths corresponding to A,C,G,T.
		    
		Returns:
		    Nothing
		    
		"""
		self.roi2psn = dict()       # region of interest -> genomic position
		self.psn2roi = dict()       # genomic position -> region of interest
		self.infotag = infotag	
		self.expectedErrorRate = expectedErrorRate
		
	def add_roi(self, roi_name, roi_positions):
		""" adds a region of interest, which is a set of genomic positions for which the
		variation should be extracted.
		
		self.roi2psn and self.psn2roi allow in memory lookup between positions and regions, and v/v;
		add_roi creates entries in roi2psn and psn2roi for all roi_names and roi_positions.
		
		Note that the roi_positions must be 1-indexed.  A value error is raised if a zero position is added.
		
		Args:
			roi_name: the name of the region of interest, example 'gene3'
			roi_position: a list containing the one indexed positions  of the bases in roi_name
			
		Returns:
			nothing
		"""
		
		try:
			self.roi2psn[roi_name] = set([])
		except KeyError:
			# already exists
			pass
		for roi_position in roi_positions:
			if roi_position == 0:
				raise ValueError("Positions supplied must be 1, not zero indexed")
			
			self.roi2psn[roi_name].add(roi_position)
			if roi_position in self.psn2roi.keys():
				self.psn2roi[roi_position].add(roi_name)
			else:
				self.psn2roi[roi_position] = set([roi_name])
				
	def persist(self, outputfile, mode='w'):
		""" persists any bases examined (i.e. which are part of rois) to an indexed hdf5 file.

		For M. tuberculosis, which has a 4.4e6 nt genome, this takes up about 52MB if
		all bases are part of ROIs (as, for example, implemented by regionScan).
		
		The bases examined are indexed by ROI and position,
		allowing near instantaneous access from on-disc stores.
		
		The HDF store access is implemented via Pandas and PyTables.
		Each matrix stored is associated with a key.
		The key used is self.guid, which is set by self.parse().
		If guid is not set, an error is raised.
		
		By default, any existing HDF file will be overwritten (mode 'w').
		To append data to an existing hdf file, use mode 'a'.
		
		Parameters:
			outputfile: the outputfile name
			mode: 'a' to append to an existing file; 'w' to overwrite.
			
		Returns:
			None
		
		"""
		
		if self.guid is None:
			raise ValueError("Cannot write hdf file with a null table key.  You must set the guid when you .parse() the vcf file")
		if self.bases is None:
			raise ValueError("No data to write; you must load data using .parse() first.")
		
		warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)		# guids are not valid python names, but this doesn't matter
		self.bases.to_hdf(outputfile,
						  key = self.guid,
						  mode=mode,
						  format='t',
						  complib='blosc:blosclz',
						  data_columns = ['roi_name', 'pos'],
						  complevel=9)
	
	def parse(self, vcffile, guid=None):
		""" parses a vcffile.
		stores a pandas data frame, with one row per roipsn/roiname combination, in self.bases.
		You must provide a guid is you wish to persist the object using the .persist method.
		Wrapper around ._parse().
		
		Arguments:
			vcffile: the vcf file to read
			guid: a guid identifier for the parsed object; required only if using .persist() to store the parsed object.
			
		Returns:
			None
		"""
		self.guid = guid
		self._parse(vcffile)
		
	def _parse(self, vcffile):
		""" parses a vcffile.
		stores a pandas data frame, with one row per roipsn/roiname combination, in self.bases
		
		Arguments:
			vcffile: the vcf file to parse
			
		Returns:
			None
			Output is stored in self.bases
		"""


		# set up variable for storing output
		resDict = {}
		nAdded = 0						
		self.region_stats = None		 
		self.bases = None
		
		# transparently handle vcf.gz files.
		if vcffile.endswith('.gz'):
			f = gzip.open(vcffile, "rb")
		else:
			f = file.open(vcffile, "rt")			

		# precompute a sorted list of positions to look for
		# in an efficient data structure
		sought_psns = deque(sorted(self.psn2roi.keys()))

		try:
			sought_now = sought_psns.popleft()

			# iterate over the vcf file
			for line in f:
				line=line.decode()
				if line[0] == "#":
					continue  # it is a comment; go to next line;
				if "INDEL" in line:
					continue  #this is not needed because ours don't contain INDELs anymore; go to next line;
				
				# parse the line.
				chrom, pos, varID, ref, alts, score, filterx, infos, fields, sampleInfo = line.strip().split()
				pos = int(pos)

				if pos == sought_now:

					alts = alts.split(",")
					infos = dict(item.split("=") for item in infos.split(";"))
					
					# confirm the self.infos tag is present.
					try:
						baseCounts4=list(map(int, infos[self.infotag].split(",")))   #get frequencies of high quality bases
					except KeyError:
						raise KeyError("Expected a tag {0} in the 'info' component of the call file, but it was not there.  Keys present are: {1}".format(self.infotag, infos.keys()))
						
					# extract the baseCounts, and do QC
					baseFreqs = baseCounts4.copy()
					baseFreqs.sort(reverse=True)  
					if not len(baseFreqs)==4:
						raise TypeError("Expected tag {0} to contain 4 depths, but {1} found.  Base = {2}; tag contents are {4}".format(self.infos,len(baseFreqs), pos, baseCounts4))
					depth = sum(baseCounts4)
					
					# compute probability that the minor variant frequency differs from self.expectedErrorRate from exact binomial test
					if (baseFreqs[0]<depth and depth>0):        # the majority base is not the only base AND depth is more than 0;
						pvalue=stats.binom_test(x=baseFreqs[1],n=depth,p=self.expectedErrorRate)   # do the test if any variation
					elif baseFreqs[0]==depth:
						   pvalue=1        # there is only one base
					elif depth==0:
						   pvalue=None     # can't tell, no data
					else:
						   raise Error("Logical error: should never reach this point {0} {1}".format(baseFreqs[0], depth))
	
					if pvalue==0:
						mlp= 250        # code minus log p as 250 if p value is recorded as 0 in float format
					elif pvalue is not None:
						mlp= -math.log(pvalue,10)
					elif pvalue is None:
						mlp=None
						   
					# store output in a dictionary   
					if depth>0:
						maf=float(baseFreqs[1])/float(depth)
					else:
						maf=None
		
					for roi_name in self.psn2roi[sought_now]:
						nAdded += 1
						resDict[nAdded] = {'roi_name':roi_name, 'pos':pos, 'ref':ref, 'depth':depth,\
									'base_a':baseCounts4[0],
									'base_c':baseCounts4[1],
									'base_g':baseCounts4[2],
									'base_t':baseCounts4[3], \
									'maf':maf,
									'mlp':mlp}                                 

					# recover the next item to recover
					try:
						sought_now = sought_psns.popleft()
					except IndexError:		# no positions selected
						break				# all positions have been selected
					
		except IndexError:		# no positions defined for selection; this is allowed
			pass
		
		# construct data frame
		self.bases=pd.DataFrame.from_dict(resDict, orient='index')
	
		# construct summary by region, defined by roi_name
		if len(self.bases.index)>0:
			r1= self.bases.groupby(['roi_name'])['depth'].mean().to_frame(name='mean_depth')
			r2= self.bases.groupby(['roi_name'])['depth'].min().to_frame(name='min_depth')
			r3= self.bases.groupby(['roi_name'])['depth'].max().to_frame(name='max_depth')
	
			r4= self.bases.groupby(['roi_name'])['pos'].min().to_frame(name='start')                
			r5= self.bases.groupby(['roi_name'])['pos'].max().to_frame(name='stop')                
			r6= self.bases.groupby(['roi_name'])['pos'].count().to_frame(name='length')
			
			# if all mafs are NA, then mean() will fail with a pandas.core.base.DataError
			try:                
					r8= self.bases.groupby(['roi_name'])['maf'].mean().to_frame(name='mean_maf')     
			except pd.core.base.DataError:
					r8=r1.copy()
					r8.columns=['mean_maf']
					r8['mean_maf']=None
	
			# compute total depth
			r9= self.bases.groupby(['roi_name'])['depth'].sum().to_frame(name='total_depth')
			
			# compute total_nonmajor_depth
			self.bases['most_common'] = self.bases[['base_a','base_c','base_g', 'base_t']].max(axis=1)
			self.bases['nonmajor'] = self.bases['depth'] - self.bases['most_common']
			r10= self.bases.groupby(['roi_name'])['nonmajor'].sum().to_frame(name='total_nonmajor_depth')
									
			df=pd.concat([r1,r2,r3,r4,r5,r6,r8, r9, r10], axis=1)              # in R,  this is a cbind operation
		else:
			df= None
		self.region_stats= df
		f.close()
		
class regionScan_from_genbank(vcfScan):
	""" annotates regions of the genome from a genbank entry.
	
	parses a vcf file, extracting high-quality bases at regions of the genome
	defined from the features of a genbank entry corresponding to the genome to which the
	data was mapped.
	
	The regions can be defined in two ways:
		* CDS - annotate CDS and intergenic regions
		* All - annotate just one region, covering the entire genome.
	
	For each region, computes minor variant frequencies for each set of lineage defining positions.
	For example, if 500 positions define the open reading from of gene B55, it will compute the minor variant frequency
	across these 500 positions, as (total minor variant frequency across the 500 bases)/ (total depth across the 500 bases)
	
	This class inherits from vcfScan, which provides much of the functionality.
	It provides new __init__ and .parse methods relative to vcfScan.
	The __init__ method (constructor) defines the regions of interest from the genbank file.
	The the .parse method computes for each region the minor allele frequency.
	
	Example usage:
	
	# read the genbank file; define regions based on CDS	
	genbank_file_name = os.path.join("..", "testdata", "NC_0103971.1.gb")
	rs = regionScan_from_genbank(genbank_file_name, method = 'CDS')
	
	# export the extracted CDs to excel;
	output_excel = os.path.join('..','testdata', 'regions.xlsx')
	rs.regions.to_excel(output_excel)

	# parse a vcf file, analysing all bases
	vcfname = os.path.join("..", "testdata", 'NC_0103971_example1.vcf.gz')
	rs.parse(vcfname)

	# optional: store the analysed genome in an indexed, compact (~ 30% of vcf.gz file size) format
	# this is a useful activity iff it might need to be analysed subsequently.
	outputfile = os.path.join('..','unitTest_tmp','test.h5')
	rs.persist(outputfile=outputfile)
	
	# export an analysis of all regions
	guid = 'test'
	region_stats_filename = os.path.join("..", "testdata", 'test.csv')
	rs.region_stats.to_csv(region_stats_filename)
	
	"""
	
	def __init__(self,
				 genbank_file_name,
				 method = 'CDS',
				 expectedErrorRate = 0.001,
				 infotag = 'BaseCounts4',
				 min_region_size = 15):
		""" defines the regions to study based on a
		genbank file.
		
		Important: internally, positions delivered are ** 1 indexed ** like the VCF file,
		not like BioPython which zero-indexes positions.
		
		Arguments:
			genbank_file_name: NCBI genbank format file used to define the regions of interest
			method: either 'CDS' or 'all'
			* method = 'CDS' - annotate CDS and intergenic regions.
							 For example, if there are genes G1, G2, and G3
							 with CDS 100..1000, 1200..2000, and 2200..3000,
							 
							 the software will define regions
							 
							 1..99 'near_G1'
							 100..1000 'G1'
							 1001..1199 'near_G2'
							 1200..2000 'G2'
							 2001..2199 'near_G3'
							 2200..3000 'G3'
			* method = 'All' - annotate just one region, covering the entire genome.
		    * expectedErrorRate: the expected per-base error.  if q30 filters are use, the expectedErrorRate is <= 0.001.
		    * infotag: the tag in the VCF INFO section containing the base counts to analyse.
		    * min_region_size: the minimum size of a region.
			                If regions identified are less than min_region_size in length, then they are combined with the next region.
		Output:
			Stores features in self.regions, as a Pandas dataframe.
			The positions delivered are ** 1 indexed ** like the VCF file, not like BioPython
	
			You can export these using
			[object].regions.to_excel(file='filename.xlsx')
			
		Returns:
			None
		"""
		self.expectedErrorRate = expectedErrorRate
		self.roi2psn = dict()
		self.psn2roi = dict()
		self.regions = None
		self.region_stats = None
		self.infotag = infotag
		# read data from genbank
		record = SeqIO.read(genbank_file_name, format='genbank')
	
		# extract features
		if method == 'All':
			self.add_roi(roi_name = 'All', roi_positions = range(1, len(record.seq)+1))
			feature_list = [{'id':1,
							 'name':'all',
							 'start_pos':1,
							 'end_pos':len(record.seq)}
				]
			
		elif method == 'CDS':
			# extract the CDs and intergenic regions.
			previous_end = -1
			feature_id = -1
			feature_list = []

			for feature in record.features:
				if feature.type == 'CDS':
	
					featname = feature.qualifiers.get('locus_tag')
					nt = set()
					for location in feature.location:		# iterates over each base
						nt.add(location)
					
					cds_name = featname[0]	
					current_start = min(nt)
					current_end = max(nt)
		
					# if there is a gap between this and the previous feature,
					# then mark as an intergenic region;
					if not current_start == previous_end +1:
						feature_id +=1
						feature_list.append({'id':feature_id, 'name':'near_'+cds_name,'start_pos':previous_end+1+1, 'end_pos':current_start-1+1})
						
					# add the current feature	
					feature_id +=1	
					feature_list.append({'id':feature_id, 'name':cds_name,'start_pos':current_start+1, 'end_pos':current_end+1})
		
					previous_end = current_end
				
			if not previous_end == len(record.seq):		# there's an extra bit at the right hand end
				feature_id +=1
				feature_list.append({'id':feature_id, 'name':'near_end', 'start_pos':previous_end+1+1,  'end_pos':len(record.seq)})
		
			# now iterate over feature_list, merging any features with length < min_region_length
			features_final = {}
			features_merged ={}
			for feat in feature_list:
				if np.abs(feat['end_pos']-feat['start_pos']) < min_region_size:
					features_merged[feat['id']] = feat
				else:
					features_final[feat['id']] = feat

			for id in features_merged.keys():
				merge_to = id
				while merge_to>0:
					merge_to = merge_to -1
					if merge_to in features_final.keys():
						features_final[merge_to]['start_pos'] = min(features_final[id-1]['start_pos'],features_merged[id]['start_pos'])
						features_final[merge_to]['end_pos'] = max(features_final[id-1]['end_pos'],features_merged[id]['end_pos'])
						if not features_final[merge_to]['name'][0:4] == 'incl':
							features_final[merge_to]['name'] = 'incl_' + features_final[merge_to]['name']
						break
	
			feature_list = []
			for i in features_final.keys():
				region = features_final[i]
				feature_list.append(region)
				self.add_roi(roi_name = region['name'], roi_positions = range(region['start_pos'], region['end_pos']+1))
		
		else:
			raise ValueError("Only methods 'CDS' and 'All' are supported; got {0}".format(method))
		
		# create an easy to display list of what we are studying
		# Creating internal region spreadsheet.  Export as .regions.to_excel(filename)")

		
		self.regions = pd.DataFrame.from_records(feature_list, index='id')
		for i in self.regions.index:
			self.regions.loc[i,'length']= np.abs(self.regions.loc[i,'end_pos']-self.regions.loc[i,'start_pos'])
			self.regions.loc[i,'mid']= (self.regions.loc[i,'end_pos']+self.regions.loc[i,'start_pos'])	/2		

		
class test_vcfScan_1(unittest.TestCase):
	def runTest(self):
		""" tests definition of regions """
		v = vcfScan()
		v.add_roi('One',set([1,2,3]))
		self.assertEqual(v.roi2psn, {'One':set([1,2,3])})
		self.assertEqual(v.psn2roi,
						 {1:{'One'}, 2:{'One'}, 3:{'One'}})

		v.add_roi('Two',set([2,3,4]))
		self.assertEqual(v.roi2psn, {'One':set([1,2,3]), 'Two':set([2,3,4])})
	
		self.assertEqual(v.psn2roi,
						 {1:set(['One']), 2:set(['One','Two']), 3:set(['One', 'Two']), 4:set(['Two'])})

		with self.assertRaises(ValueError):
			v.add_roi('Not allowed', set([0]))
		
class test_vcfScan_2(unittest.TestCase):
	def runTest(self):
		""" tests reading from a region when none is specified """
		v = vcfScan()
		v.add_roi('One',set([]))
		
		inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
		
		if not os.path.exists(inputfile):
			self.fail("Input file does not exist.  Please see README.  You may need to install test data.")
		v.parse(vcffile = inputfile)
		self.assertEqual(len(v.bases.index),0)

class test_vcfScan_3(unittest.TestCase):
	def runTest(self):
		""" tests reading when the info tag does not exist"""

		v = vcfScan(infotag = 'missing')
		inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
		if not os.path.exists(inputfile):
			self.fail("Input file does not exist.  Please see README.  You may need to install test data.")
			
		# should raise an error
		v.add_roi('One', [1,2,3])
		
		with self.assertRaises(KeyError):
			v.parse(vcffile = inputfile)

class test_vcfScan_4(unittest.TestCase):
	def runTest(self):
		""" tests reading from a region """
		v = vcfScan()
		v.add_roi('One',set([1,2,3]))
		v.add_roi('Two',set([2,3,4]))
		inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
		if not os.path.exists(inputfile):
			self.fail("Input file does not exist.  Please see README.  You may need to install test data.")
		v.parse(vcffile = inputfile)
		self.assertEqual(len(v.bases.index),6)

class test_regionScan_1a(unittest.TestCase):
	def runTest(self):
		
		genbank_file_name = os.path.join("..", "testdata", "NC_0103971.1.gb")
		if not os.path.exists(genbank_file_name):
			self.fail("Input file does not exist.  Please see README.  You may need to install test data.")
		rs = regionScan_from_genbank(genbank_file_name, method= 'CDS', min_region_size = 0)

		self.assertTrue(isinstance(rs.regions, pd.DataFrame))
		self.assertEqual(len(rs.regions.index), 9795)
		
		# check nothing is missing from our regions, and nothing is covered twice
		bases_covered = set()
		for i in rs.regions.index:
			nt = range(rs.regions.loc[i,'start_pos'], rs.regions.loc[i,'end_pos']+1)
			for pos in nt:
				#if pos in bases_covered:
				#	self.fail("Duplicate base (overlapping regions) {0}".format(pos))		#overlapping cds actually exist, e.g. 4275
				bases_covered.add(pos)
		seqlen = 5067172
		bases_expected = set(range(1,seqlen+1))
		self.assertEqual(bases_expected, bases_covered)

		n = len(rs.regions.query('length<15').index)
		self.assertTrue(n>0)

class test_regionScan_1b(unittest.TestCase):
	def runTest(self):
		
		genbank_file_name = os.path.join("..", "testdata", "NC_0103971.1.gb")
		if not os.path.exists(genbank_file_name):
			self.fail("Input file does not exist.  Please see README.  You may need to install test data.")
		rs = regionScan_from_genbank(genbank_file_name, method= 'CDS', min_region_size=15)
		
		self.assertTrue(isinstance(rs.regions, pd.DataFrame))
		self.assertTrue(len(rs.regions.index)< 9795)
		
		# check nothing is missing from our regions, and nothing is covered twice
		bases_covered = set()
		for i in rs.regions.index:
			nt = range(rs.regions.loc[i,'start_pos'], rs.regions.loc[i,'end_pos']+1)
			for pos in nt:
				#if pos in bases_covered:
				#	self.fail("Duplicate base (overlapping regions) {0}".format(pos))		#overlapping cds actually exist, e.g. 4275
				bases_covered.add(pos)
		seqlen = 5067172
		bases_expected = set(range(1,seqlen+1))
		self.assertEqual(bases_expected, bases_covered)
		n = len(rs.regions.query('length<15').index)
		self.assertTrue(n == 0)
		print(rs.regions)
		
class test_regionScan_2(unittest.TestCase):
	def runTest(self):
		
		genbank_file_name = os.path.join("..", "testdata", "NC_0103971.1.gb")
		if not os.path.exists(genbank_file_name):
			self.fail("Input file does not exist.  Please see README.  You may need to install test data.")
			
		rs = regionScan_from_genbank(genbank_file_name, method = 'All')
		
		self.assertTrue(isinstance(rs.regions, pd.DataFrame))
		self.assertEqual(len(rs.regions.index), 1)
		
		# check nothing is missing from our regions, and nothing is covered twice
		bases_covered = set()
		for i in rs.regions.index:
			nt = range(rs.regions.loc[i,'start_pos'], rs.regions.loc[i,'end_pos']+1)
			for pos in nt:
				#if pos in bases_covered:
				#	self.fail("Duplicate base (overlapping regions) {0}".format(pos))		#overlapping cds actually exist, e.g. 4275
				bases_covered.add(pos)
		seqlen = 5067172
		bases_expected = set(range(1,seqlen+1))
		#print(max(bases_expected))		# genome length is 5067172
		self.assertEqual(bases_expected, bases_covered)

			
class test_regionScan_3(unittest.TestCase):
	def runTest(self):
		""" tests persistence of the bases analysed to hdf5 """
		genbank_file_name = os.path.join("..", "testdata", "NC_0103971.1.gb")
		rs = regionScan_from_genbank(genbank_file_name, method = 'All')
			
		guid = 'test_guid'
		inputfile=os.path.join("..",'testdata','NC_0103971_example1.vcf.gz')
		if not os.path.exists(inputfile):
			self.fail("Input file does not exist.  Please see README.  You may need to install test data.")
		res = rs.parse(vcffile = inputfile, guid = guid)
	
		outputfile = os.path.join('..','unitTest_tmp','test.h5')
		if os.path.exists(outputfile):
			os.unlink(outputfile)		# delete it
		rs.persist(outputfile=outputfile)
	
		df = pd.read_hdf(outputfile, guid, where = ['pos in [2,3,6,1500,15000,150000, 1000000, 5067172]'])
		self.assertEqual(len(df.index),8)
		df = pd.read_hdf(outputfile, guid, where = ['roi_name=="MAB_0992"'])
		self.assertEqual(len(df.index), 0)

		seqlen = 5067172
		df = pd.read_hdf(outputfile, guid, where = ['roi_name=="All"'])
		self.assertEqual(len(df.index), seqlen)
		
		
		# now test the results of the persistence operation;
		df = pd.read_hdf(outputfile, guid)

		# check every base has been examined. 
	
		bases_expected = set(range(1,seqlen+1))
		bases_observed = set(df['pos'])
		self.assertEqual(len(bases_observed - bases_expected),0)
		self.assertEqual(len(bases_expected - bases_observed),0)

