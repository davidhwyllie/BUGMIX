#!/usr/bin/python

# pyssandra 
import csv
import os
import sys
import logging
import shutil
import glob
import unittest
import re
from Bio import SeqIO	
import gzip
from pyssandra import *
import readline
import time
import pyssandra
import uuid
import datetime

logging.getLogger().setLevel(logging.WARNING)

class fastaComposition():
	""" returns the base composition of a fasta file
	Both gzipped and uncompressed fasta files are supported.
	
	usage:
	fc=fastaCompostion()
	for composition in fc.composition(filename):
		print composition
		"""
	def __init__(self):
		pass	# does nothing
	def composition(self,filename):
		if filename.endswith('.gz'):
			handle = gzip.open(filename)
		else:
			handle=open(filename)
		with handle:
			for seq_record in SeqIO.parse(handle, "fasta"):
				
				retDict=dict()
				retDict['length']=len(seq_record)
				for base in ('A','C','G','T','N'):
					retDict[base+'_count']=seq_record.seq.count(base)
					if retDict['length']>0:
						retDict[base+'_prop']=float(retDict[base+'_count'])/float(retDict['length'])
		
				yield retDict
class Test_fastaComposition_1(unittest.TestCase):
	""" tests input validation """
	def runTest(self):
		fc=fastaComposition()
		nreturned=0
		for composition in fc.composition('test.fasta.gz'):
			nreturned=nreturned+1

		self.assertEqual(nreturned,1)
		nCount=composition['N_count']
		self.assertEqual(nCount,349509)	
class Test_fastaComposition_1(unittest.TestCase):
	""" tests input validation """
	def runTest(self):
		fc=fastaComposition()
		nreturned=0
		for composition in fc.composition('../testdata/test2.fasta'):
			nreturned=nreturned+1			
		self.assertEqual(nreturned,5)
class CassandraDataReader():
	""" methods for reading data from Cassandra
	
	cdr=CassandraDataReader()
	cdr.addSamplesByReference()
	for sample in cdr.samples():
		print "SampleId=",str(sample)		# the guid
		print cdr.sampleAnnotation(sample)
		print cdr.samplePath(sample)
	
	# thanks to Oriol for parts of this code
	"""
	def __init__(self, fileRoot='/mnt/microbio/ndm-hicf/ogre'):
		""" constructor """
		self.refresh()
		self.q=CsQuery()
		self.fileRoot=fileRoot
		
	def _annotationDictionary(self,obj):
		""" if passed a CsElement or other compatible object, returns attributes as a dictionary """
		# no unit tests
		retDict={}
		for attribute_name in obj._entity.attr_names:
			content=getattr(obj, attribute_name)

			# it's perhaps a dictionary
			if type(content) is dict:
				for key in content.keys:
					composite_key=attribute_name+"_"+key
					if not content[key] is None:
						retDict[composite_key]=content[key]
						
			else:	# coerce to a string
				if not content is None:
				   retDict[attribute_name]=content   
		return(retDict)
	def plates(self):
		for thisplate in self._plates:
			yield thisplate		
	def studies(self):
		for thisstudy in self._studies:
			yield thisstudy			
	def samples(self):
		for thissample in self._samples:
			yield thissample
	def samplePath(self, sample_id):
		""" returns a path to a file type, in this case fasta """
		pathToFile= "{0}/pipeline_output/{1}/MAPPING/{2}_{3}/STD/basecalls/{4}_v3.fasta.gz".format(self.fileRoot, sample_id,self._reference_id,self._reference_name, sample_id)
		if self._reference_id is None or self._reference_name is None:
			raise NotImplementedError		# self._reference_id or self._reference_name is none;
		else:
			return pathToFile
	def sampleAnnotation(self, sample):
		""" pulls descriptive data about a sample out into a list of dictionaries """
		# no unit tests
		annotations={}

		sampleObj=self.q.get('sample',filterop="id={0}".format(str(sample))) # get the sample object
		sampleObj=sampleObj.fetch()[0]
		annotations['csSample']=self._annotationDictionary(sampleObj)
		
		plates=sampleObj.is_related_with_plate
		plates=plates.fetch()
		for plate in plates:
			annotations['csPlate']=self._annotationDictionary(plate)
	
		# add the study
		studies=sampleObj.is_related_with_study
		studies=studies.fetch()
		for study in studies:
			annotations['study']=self._annotationDictionary(study)
			
		# recover mapcall information if a reference has been specified
		# get the mapcall objects; there may be multiple references, so we just select the one we are using		
		for id in sampleObj.has_mapcall:	
			mapcallObj=self.q.get('mapcall',filterop="id={0}".format(id))			
			mapcallObj=mapcallObj.fetch()[0]
			annotDict=self._annotationDictionary(mapcallObj)
			if 'contig_id' in annotDict:
				if annotDict['contig_id']==self._reference_name:
					annotations['mapcall']=annotDict
		
		x=sampleObj.is_related_with_study.fetch()[0]		# this adds one and only one study
		annotations['study']=self._annotationDictionary(x)	
					
		return(annotations)		
	def refresh(self):
		""" removes any data selected"""
		self._studies=set()
		self._plates=set()
		self._samples=set()
		self._reference_name=None
		self._reference_id=None
		
	def addAllStudies(self):
		""" adds all studies to the object"""
		res=self.q.get('study')
		for i in res:
			self._studies.add(i)
	def addStudiesByStudyCode(self, studyCd):
		""" adds any studies including study code, e.g. Blood-culture-extractions"""
		res=self.q.get('study', filterop='name=%s' % studyCd)
		for i in res:
			self._studies.add(i)
	def addSamplesByReference(self, reference_name="R00000039"):
		""" gets all sampleIds mapped to a reference """

		### get reference
		reference = self.q.get('reference',filterop="name={0}".format(reference_name))
		self._reference_name=reference_name
		self._reference_id= reference.fetch()[0].id
		
		### new code: use two search strategies as these have been used by slightly different pipeline versions in the past
		# strategy 1
		samples = reference.is_related_with_organism.is_related_with_isolate.is_related_with_sample
		for sample in samples:
			self._samples.add(sample)

		###get samples, strategy 2
		samples = reference.is_related_with_organism.is_related_with_speciation.is_related_with_sample
		for sample in samples:
			self._samples.add(sample)
				
class Test_CassandraDataReader_1(unittest.TestCase):
 	""" tests pyssandra: whether any studies can be found """
 	def runTest(self):
		sm=CassandraDataReader()
		sm.addAllStudies()
		nStudies=0
		for i in sm.studies():
			nStudies=nStudies+1
		self.assertEqual(nStudies>0, True)
class Test_CassandraDataReader_2(unittest.TestCase):
 	""" tests pyssandra: whether any studies can be found """
 	def runTest(self):
		sm=CassandraDataReader()
		sm.addStudiesByStudyCode('MtubPilot')		# should exist
		nStudies=0
		for i in sm.studies():
			nStudies=nStudies+1
		self.assertEqual(nStudies, 1)		
class Test_CassandraDataReader_3(unittest.TestCase):
 	""" tests pyssandra: whether any studies can be found """
 	def runTest(self):
		sm=CassandraDataReader()
		sm.addStudiesByStudyCode('zxcvzxcvzxcsdgdgfg')		# should not exist
		nStudies=0
		for i in sm.studies():
			nStudies=nStudies+1
		self.assertEqual(nStudies, 0)
class Test_CassandraDataReader_4(unittest.TestCase):
 	""" tests pyssandra: whether any sets are empty when created """
 	def runTest(self):
		sm=CassandraDataReader()
		n=len(sm._plates)
		self.assertEqual(n, 0)
class Test_CassandraDataReader_5(unittest.TestCase):
 	""" tests pyssandra: whether any sets are empty when created """
 	def runTest(self):
		sm=CassandraDataReader()
		n=len(sm._samples)
		self.assertEqual(n, 0)
class Test_CassandraDataReader_6(unittest.TestCase):
 	""" tests pyssandra: whether any reads mapped to TB are present """
 	def runTest(self):
		sm=CassandraDataReader()
		sm.addSamplesByReference(reference_name="R00000039")
		n=len(sm._samples)
		self.assertEqual(n>0, True)
class Test_CassandraDataReader_7(unittest.TestCase):
 	""" tests pyssandra: whether any sets are empty when created """
 	def runTest(self):
		sm=CassandraDataReader()
		n=len(sm._plates)
		self.assertEqual(n==0, True)

# run unittests
if __name__ == '__main__':
	logging.getLogger().setLevel(logging.WARNING)
	unittest.main()            ## test everything
