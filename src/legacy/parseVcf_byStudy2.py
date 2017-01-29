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
from scipy import stats
from numpy import median, array, mean, std
import time
from csv import DictWriter
import math
import cPickle
from CassandraAndOgreAccess import *
logging.getLogger().setLevel(logging.CRITICAL)

class vcfParser():
    """ methods for extracting base counts from vcf file
    
    optionally excludes any bases in the maskedBases() set passed to the constructor
    
    example usage:
    
    print time.strftime("%c")            
            inputfile='/home/dwyllie/TBMIX/real/testvcf/BG134_Mtub/MAPPING/R00000039/STD/basecalls/5229afb0-879f-4989-b541-262de81f7483_R00000039_v3.vcf.gz'
    vp=vcfParser(inputfile, '5229afb0-879f-4989-b541-262de81f7483')
    
    # optionally can change a series of defaults
    vp.includeAll=True    # output for each base, much slower
    vp.maxDepthDistribution=200   # depth histogram only up to 200
    vp.maxLines=1e3       # only the first 1000 for testing.
    
    vp.evaluateBaseCalls()  # do the main calculation
    
    # now output 
    print vp.depthDistribution
    print vp.depthDistributionDictionary()
    print vp.unclusteredPositions()
    print vp.minorAlleleFrequencies()
    print vp.minorAlleleDictionary()
    print vp.summaryStatistics()
    print time.strftime("%c")

    """
    def __init__(self, vcffile, sampleId, badBases=set()):
        """ constructor """
        self.filename=vcffile
        self.sampleId=sampleId          # goes into the output
       
        ## parameters
        self.includeAll=False   	# do not report invariant sites if False
        self.maxDepthDistribution=1000   # compute depth histograms up to 1000 deep
        self.expectedErrorRate=0.001    # corresponds to q30
        self.qcutoff=0.01               # after fdr, if q<0.01 then regard this as significant
        self.maxLines=1e10              # 1e10 for everything
        self.badBases=badBases		# masked bases which are not reported
	
        ## define data structures used during the parse operation
        self.timings={}
        self.depthDistribution=[0]*self.maxDepthDistribution
        self.unadjustedpValues={}
        self.baseproperties={}
        self.significantPositions=[]   # list to hold sites of significant variation
        self.clusteredPositions=[]     # list to hold sites which are clustered together

    def exportDictionary(self, DictofDicts, targetfile, format='csv'):
        """ writes one of the internal dictionary of dictionaries to file, in a variety of formats.
        Internal function.  At present, only csv is supported.
        The dictionaries are written out in the order of the keys.
        The dictionaries which represent each row are written out using csv.DictWriter.
        """
        
        # create the directory if it doesn't exist
        self.ensuredir(os.path.dirname(targetfile))
        # scan all dictionaries and obtain all the keys within them.
        fields=set()
        for thisKey in sorted(DictofDicts):
            for thisSubKey in DictofDicts[thisKey]:
                if thisSubKey not in fields:
                    fields.add(thisSubKey)
        fieldnames=sorted(fields)
        with open(targetfile,'w') as csvfile:
            writer=DictWriter(csvfile,fieldnames=fieldnames)
            writer.writeheader()
            for thisKey in sorted(DictofDicts):
                writer.writerow(DictofDicts[thisKey])          
    def ensuredir(self, newpath):
        """ makes a directory if it does not exist """
        if not os.path.exists(newpath):
            os.makedirs(newpath)
    def exportAll(self,basedir):
        """ writes depth distribution, summary and variant alleles to disk, and pickles the object"""
        targetdir=os.path.join(basedir,self.sampleId)
        self.exportDictionary(self.depthDistributionDictionary(),os.path.join(targetdir,'depthDistribution.csv'))
        self.exportDictionary(self.minorAlleleDictionary(), os.path.join(targetdir,'variantSites.csv'))
        self.exportDictionary(self.summaryStatistics(),os.path.join(targetdir,'summary.csv'))
        with open(os.path.join(targetdir,'pickled'),'wb') as picklefile:
           cPickle.Pickler(picklefile,2).dump(self)
    def summaryStatistics(self):
        """ outputs mean, sd, number and median of variant sites """
        x=array(self.minorAlleleFrequencies())
        return {'summary':{'sampleId':self.sampleId,'median':median(x),'mean':mean(x),'std':std(x), 'nvariants':len(x)}}    
    def unclusteredPositions(self):
        """ returns list of positions of mixed bases which are not clustered"""
        x=list(set(self.significantPositions)-set(self.clusteredPositions))
        x.sort()
        return x    
    def minorAlleleFrequencies(self):
        """ returns a vector of minor allele frequencies, as floating point numbers """
        maf=[]
        for pos in self.unclusteredPositions():           
            maf.append(float(self.baseproperties[pos]['freq_2'])/float(self.baseproperties[pos]['depth']))
        return maf    
    def minorAlleleDictionary(self):
        """ returns a dictionary of bases containing unclustered minor alleles.
        The key is the position, and the dictionary contains sampleId, read counts, and maf.
        Suitable for export using DictWriter"""    
        maf={}
        for pos in self.unclusteredPositions():
		if self.baseproperties[pos]['qvalue']==0:
			mlp= 250	# if q value recorded as 0, register minus log p as 250
		else:
			mlp= -math.log(self.baseproperties[pos]['qvalue'],10)
			maf[pos]={
				'sampleId':self.sampleId, \
				'base':pos, \
				'ma_reads':self.baseproperties[pos]['freq_2'], \
				'depth':self.baseproperties[pos]['depth'], \
				'maf':float(self.baseproperties[pos]['freq_2'])/float(self.baseproperties[pos]['depth']),
				'minusLogQvalue':mlp,
				'base_a':self.baseproperties[pos]['base_a'],
				'base_c':self.baseproperties[pos]['base_c'],
				'base_g':self.baseproperties[pos]['base_g'],
				'base_t':self.baseproperties[pos]['base_t']}
        return maf   
    def adjustedpValues(self,x):
        """
        Assumes a list or numpy array x which contains p-values for multiple tests
        Copied from p.adjust function from R.
        Implements Benjamini Hochberg (FDR) computation 
        """
       
        o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
        ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
        q = sum([1.0/i for i in xrange(1,len(x)+1)])
        l = [q*len(x)/i*x[j] for i,j in zip(reversed(xrange(1,len(x)+1)),o)]
        l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
        return l

    def depthDistributionDictionary(self):
        """ returns the depth distribution as a dictionary """
        dd={}
        for i in range(0,len(self.depthDistribution)-1):
            dd[i]={'sampleId':self.sampleId, 'depth':i, 'nPositions':self.depthDistribution[i]}
        return dd
    def evaluateBaseCalls(self, eliminateNeighbours=0):
        """ processes the base calls from a VCF file.
        Assumes that the 'baseCounts4' element is present, and contains >Q30 bases.
        From this, it computes
            a depth distribution
            a list of bases which have a minor allele frequency in excess of expected based on Q30 (by default)
            it also marks these, eliminating any within 12 bases of each other,
            a strategy which has been used to remove erroneous snps due to mapping errors or indels.
        """
        
        self.timings['start']=time.strftime("%c")       # monitor time taken
        
        f = gzip.open(self.filename, "r")
        returnedLines=0
        
        # iterate over the vcf file
        for line in f:
            if line[0] == "#": continue  # it is a comment
            if "INDEL" in line: continue #this is not needed because ours don't contain INDELs anymore
            
            chrom, pos, varID, ref, alts, score, filterx, infos, fields, sampleInfo = line.strip().split()
            pos = int(pos)
	    
	    if not pos in self.badBases:		# ignore masked bases
		alts = alts.split(",")
		infos = dict(item.split("=") for item in infos.split(";"))
		baseCounts4=map(int, infos['BaseCounts4'].split(","))   #get frequencies of high quality bases
		baseFreqs=map(int, infos['BaseCounts4'].split(","))     #get frequencies of high quality bases
		baseFreqs.sort(reverse=True)                            #get frequencies of high quality bases, sorted
		depth = sum(baseCounts4)
		
		if depth<(self.maxDepthDistribution-1):
		    self.depthDistribution[depth]=self.depthDistribution[depth]+1
		else:
		    self.depthDistribution[(self.maxDepthDistribution-1)]=self.depthDistribution[(self.maxDepthDistribution-1)]+1
		 
		# store output in a dictionary   
		retDict={'sampleId':self.sampleId, 'pos':pos, 'ref':ref, 'depth':depth,\
			   'base_a':baseCounts4[0], 'base_c':baseCounts4[1], 'base_g':baseCounts4[2], 'base_t':baseCounts4[3], \
			   'freq_1':baseFreqs[0], 'freq_2':baseFreqs[1], 'freq_3':baseFreqs[2], 'freq_4':baseFreqs[3]
		}
		
		# compute probability from exact binomial test
		if (baseFreqs[0]<depth and depth>0):
		    pvalue=stats.binom_test(x=baseFreqs[1],n=depth,p=self.expectedErrorRate)   # do the test if any variation
		else:
		    pvalue=1        # no point, just assign 1         
		self.unadjustedpValues[retDict['pos']]=pvalue   # store this output in a dictionary, as we have to adjust it             
		  
		# ignore if all calls are the same; otherwise write the output
		if (baseFreqs[0]<depth and depth>0) or self.includeAll==True:
		    returnedLines=returnedLines+1
		    if (returnedLines>self.maxLines):
                break       # debug setting
		    self.baseproperties[retDict['pos']]=retDict
           
        ## apply fdr     
        positions=self.baseproperties.keys()        # which positions we are analysing
        pvalues=[]                                  # extract the p values into a vector
        for position in positions:                  # for all the positions analysed             
            pvalues.append(self.unadjustedpValues[position])    # add the unadjusted value to a list
        adjustedpvalues=self.adjustedpValues(pvalues)   # and apply fdr
        
        npos=-1                                     # now write the fdr derived q values back into the hash dictionary
        for position in positions:                  # write the output back into the hash dictionary
            npos=npos+1
            self.baseproperties[position]['pvalue']=self.unadjustedpValues[position]
            self.baseproperties[position]['qvalue']=adjustedpvalues[npos]
            isSignificant=0
            if (self.baseproperties[position]['qvalue']<self.qcutoff):        # base is mixed
                isSignificant=1
                self.significantPositions.append(position)                    # keep a track of positions which are significant
            self.baseproperties[position]['isSignificant']=isSignificant      # store significant positions
            self.baseproperties[position]['isClustered']=0                    # if it is next to other variant positions, we'll set this later
        
        # determine whether the variant is clustered within eliminateNeighbours of any other
        self.significantPositions.sort()
        for i in range(0,len(self.significantPositions)-2,1):
            if (i>0):
                priorPos=self.significantPositions[i-1]
            thisPos=self.significantPositions[i]
            nextPos=self.significantPositions[i+1]
            if  ((i>0) and (thisPos-priorPos<eliminateNeighbours) or (nextPos-thisPos<eliminateNeighbours)) or (i==0 and nextPos-thisPos<eliminateNeighbours):
                self.baseproperties[thisPos]['isClustered']=1
                self.clusteredPositions.append(thisPos)
        
        
        self.timings['finished']=time.strftime("%c")       # monitor time taken
   
     
### objective is to determine the frequency of mixtures across samples
## startup

# check whether we have been passed a single argument
# if it's a number, n between 1 and 6 we examine the the nth plate.
# allows running multiple threads

print "Arguments:",sys.argv
analyseNth=1
analyseStep=1
if (len(sys.argv)==2):
	retchar=str(sys.argv[1])
	if (retchar in ('0','1','2','3')):
		analyseNth=int(retchar)
		analyseStep=4

print("Analysing every %i of %i plates" % (analyseNth,analyseStep))
	
basedir='/mnt/nfs/home/dwyllie/TBMIX-real/real_allbases/'
forceXMLReload=True

# we will ignore masked bases at this point.  We can mask them out later.
badBases=set()

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
targetdir=os.path.join(basedir, 'mixedbases')
fetcher1.setTargetDir(vcfdir)
maxPlates=1e6


nPlatesTotal=0
for thisplate in sm.plates():
	
	nPlatesTotal+=1
	if (nPlatesTotal % analyseStep)==(analyseNth-1):
		print "Analysing # ",nPlatesTotal,thisplate
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
		
		if nPlates>maxPlates:
			print("Terminated as maximal processed plates reached = ", nPlates)
			break
		
		for thisSampleId in pm.samples:
			vcffile=os.path.join(vcfdir, thisSampleId+'.vcf.gz')
			targetfile=os.path.join(targetdir,thisSampleId,'summary.csv')
			
			# if the file containing reads already exists, we don't need to do anything
			if not os.path.exists(targetfile):
				fetcher1.refresh()
				fetcher1.addGuuid(thisSampleId)
				print("Will examine Guuid = %s " % thisSampleId)
				fetcher1.getVcfs(download=True, maxDownloads=1)
				
				# check whether it found the vcf file using the fetcher
				if os.path.exists(vcffile):		
					
					print("Parsing: %s - %s" % (thisplate,thisSampleId))					
					# parse the vcf file
					vp=vcfParser(vcffile, thisSampleId, badBases)
					vp.evaluateBaseCalls()
					vp.exportAll(os.path.join(basedir,'mixedbases'))
					os.remove(vcffile)
				else:
					print("Vcf file %s does not exist" % vcffile)
	else:
		print "Ignored # ",nPlatesTotal,thisplate
					
print('Complete')
exit(0)
