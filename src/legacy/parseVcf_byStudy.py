#!/usr/bin/python

# tests whether TB sequences are mixed

import argparse
from pyOgre import OgreGw,OgreEntity
from WalrusBoto import WalrusIface
import csv
import os
import sys
import gzip
import logging
import shutil
import glob
import pysam
from lxml import etree as ET
import urllib	# response library is better but is not installed
import re
from CassandraAndOgreAccess import *

logging.getLogger().setLevel(logging.CRITICAL)

class vcfParser():
	""" methods for extracting base counts from vcf file """
	def __init__(self, vcffile):
		""" constructor """
		self.filename=vcffile
		self.includeAll=False 	# don't report invariant sites

	def basecalls(self):
		f = gzip.open(self.filename, "r")

		for line in f:
		    if line[0] == "#": continue
		    if "INDEL" in line: continue #this is not needed because ours don't contain INDELs anymore
		    
		    chrom, pos, varID, ref, alts, score, filterx, infos, fields, sampleInfo = line.strip().split()
		    pos = int(pos)
		    alts = alts.split(",")
		    infos = dict(item.split("=") for item in infos.split(";"))
		    baseCounts4=map(int, infos['BaseCounts4'].split(",")) #get frequencies of high quality bases
		    depth = sum(baseCounts4)
		    
		    # ignore if all calls are the same; otherwise write the output
		    if max(baseCounts4)<depth or self.includeAll==True:
			yield([pos, ref, depth, baseCounts4[0], baseCounts4[1], baseCounts4[2], baseCounts4[3]])

### objective is to prototype a QC application for the mtubpilot study (TB in B'ham and other parts of PHE)
## if directories aren't cleaned out, will cache bam files and plate xml since to reduce network connectivity

basedir='/mnt/nfs/home/dwyllie/TBMIX-real/real/'
forceXMLReload=False

subdirs={}
for thisdirstem in ['xmlcache','vcfcache']:
	thisdir=os.path.join(basedir, thisdirstem)
	subdirs[thisdirstem]=thisdir
	print thisdir
	try:
		os.mkdir(thisdir)
	except:
		# remove all files in the directory
		# to be safe(r), we'll only delete bam xml html and csv files
		for filetype in ['.bam','.xml','.html','.csv','.vcf.gz']:
			filelist=[f for f in os.listdir(thisdir) if f.endswith(filetype)]
			for f in filelist:
				if forceXMLReload==True and filetype=='.xml':
					os.remove(os.path.join(thisdir, f))
				elif not filetype=='.xml':
					os.remove(os.path.join(thisdir, f))

## create a cassandra loader object.
# the cassandra url is hard coded in the constructor, but can be overridden by replacing 'None'
lc=loadFromCassandra(None)

## define the correct stud(ies) to be examined
studyxmlfile=os.path.join(subdirs['xmlcache'],'studies.xml')	# define a file name to store the study xml as
lc.download('studies',[],studyxmlfile)			# download it
sm=StudyManager()					# create a study manager object 
sm.loadXMLfromFile(studyxmlfile)			# parse the study xml
sm.addStudiesByStudyCode('MtubPilot')			# select one study (if add additional lines, will add extra studies toto list)

# show plates in selected studies
pm=PlateManager()	# create plate manager
nPlates=0 		# number of plates scanned

fetcher1=WalrusFetcherByGuuid('ogre')
vcfdir=os.path.join(basedir,'vcfcache')
fetcher1.setTargetDir(vcfdir)

for thisplate in sm.plates():
	
		# download the complete xml for this plate
		targetfile=os.path.join(subdirs['xmlcache'], "%s.xml" % thisplate)
		# if this file does not already exists in our cache, or if it has zero size, then download it.
		if not os.path.exists(targetfile):
			# download it
			lc.download('fullplate',[thisplate],targetfile)
		else:
			# if it's empty, delete it
			if os.path.getsize(targetfile)==0:
				os.remove(targetfile)
				# download it
				lc.download('fullplate',[thisplate],targetfile)
	
		# parse it
		nPlates+=1
		pm.refresh()							# start a new plate
		pm.loadXMLfromFile(targetfile)					# load xml
		pm.parsePlate()							# parse it
		
		for thisSampleId in pm.samples:
			vcffile=os.path.join(vcfdir, thisSampleId+'.vcf.gz')
			targetfile=os.path.join(vcfdir,thisSampleId+'.txt')
			
			# if the file containing reads already exists, we don't need to do anything
			if not os.path.exists(targetfile):
				fetcher1.refresh()
				fetcher1.addGuuid(thisSampleId)
				fetcher1.getVcfs(download=True, maxDownloads=1)
				
				# check whether it found the vcf file using the fetcher
				if os.path.exists(vcffile):		
					
					print("%s - %s" % (thisplate,thisSampleId))
					
					# parse the vcf file
					parser=vcfParser(vcffile)	# create a parser for vcffiles
					parser.includeAll=True		# force it to output base counts for all bases, not just variants
					with open(targetfile,"w") as f1:
						csvwriter=csv.writer(f1, delimiter=',')
						for l in parser.basecalls():
							csvwriter.writerow(l)
							
					os.remove(vcffile)           

				  
print('Complete')
exit(0)
