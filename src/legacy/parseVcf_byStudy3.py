#!/usr/bin/python

# tests whether TB sequences are mixed
import csv
import os
import sys
import gzip
import shutil
import glob
import pysam
from scipy import stats
import numpy  as np
#import median, array, mean, std, int as np.int
import time
from csv import DictWriter
import math
import cPickle
import unittest
import logging
import pandas as pd
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, Float, DateTime, Boolean, MetaData, select, func, LargeBinary
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref
from sqlalchemy import exc  # exceptions


### beginning of vcfSQL definitions
db=declarative_base() # classes mapping to persistent database inherit from this; global entity
class vcfBases(db):
        """ stores depths for each base at each positions """
        __tablename__ = 'vcfBases'
        pos=Column(Integer, primary_key=True)
        ref=Column(String(1))
        depth=Column(Integer, index=True)
        base_a=Column(Integer)
        base_c=Column(Integer)
        base_g=Column(Integer)
        base_t=Column(Integer)
        maf=Column(Float, nullable=True, index=True)
        mlp=Column(Float, nullable=True)
        mlq=Column(Float, nullable=True, index=True)

class vcfSQL():
        """ manages a database backend (default: sqlite) for use for VCF storage.  
        You would not normally use this class directly.  Use vcfStore instead, which is a wrapper around vcfSQL.
        """
        def __init__(self, sampleId, persistenceDir, db=db ):
            """ creates a connection to an SQLite database called '{sampleId}.db' in persistenceDir.
            If the database does not exist, it is created.          
            Does not trap for errors if the filename, directory etc are invalid
            """
            self.connString='sqlite:///{0}/{1}.db'.format(persistenceDir,sampleId)
            self.persistenceDir=persistenceDir
            self.Base = db             
            self.engine = create_engine(self.connString)    # connection string for a database
            self.Base.metadata.create_all(self.engine)  # create the table(s) if do note exist
            Session = sessionmaker(bind=self.engine)    # class
            self.session=Session()
            
        def populate(self, df):
            """ loads the data in pandas dataframe df into the table vcfBases.
            
            df must exist and be the right format, e.g. as produce by vcfStore.
            required columns are these:
            {'pos','ref','depth','base_a','base_c','base_g','base_t','maf','mlp','mlq'}
            Others will not be written """
            current_columns=set(df.columns)
            required_columns={'pos','ref','depth','base_a','base_c','base_g','base_t','maf','mlp','mlq'}
            extra_columns=current_columns - required_columns
            df.drop(extra_columns, axis=1, inplace=True)
            df.to_sql('vcfBases',index=False, if_exists='append', con=self.engine)
            
               
class vcfStore():
    """ methods for allowing rapid access to aspects of information within a vcf file.
    
    * extracting base counts from vcf file;
    * testing whether these are compatible with a prior for sequencing error rate;
    * storing them in an indexed, rapid-access manner
    * allowing access to, and computations on, single samples within this store.
    
    ## TODO
    example usage:
    
    print time.strftime("%c")            
            inputfile='/home/dwyllie/TBMIX/real/testvcf/BG134_Mtub/MAPPING/R00000039/STD/basecalls/5229afb0-879f-4989-b541-262de81f7483_R00000039_v3.vcf.gz'
    vp=vcfStore(inputfile, '5229afb0-879f-4989-b541-262de81f7483')
    
    # optionally can change a series of defaults
    vp.includeAll=True    # output for each base, much slower
    vp.maxDepthDistribution=200   # depth histogram only up to 200
    vp.maxLines=1e3       # only the first 1000 for testing.
    
    vp.evaluateBaseCalls()  # do the main calculation
    
    # now output 
    print vp.depthDistribution
    print vp.depthDistributionDictionary()
    print vp.minorAlleleFrequencies()
    print vp.alleleDictionary()
    print vp.summaryStatistics()
    print time.strftime("%c")

    """
    def __init__(self, persistenceDir, maxLines=1e8, maxDepthDistribution=1000, expectedErrorRate=0.001):
        """ includeAll                : if false, do not report invariant sites in the output
            maxDepthDistribution      : compute depth histograms up to this deep
            expectedErrorRate         : the expected error rate from sequencing; used for negative binomial test; 0.01 corresponds to q30
            qcutoff                   : labels 'significant' if FDR corrected significance is less than or equal to thisdefault: 0.01, which after fdr, if q<0.01 then regard this as significant
            maxLines                  : lines to read.  Set to 1000 for unittesting; set to 1e10 for everything
       """
        
        ## parameters
        self.maxDepthDistribution=maxDepthDistribution      # compute depth histograms up to 1000 deep
        self.expectedErrorRate=expectedErrorRate            # 0.01 corresponds to q30
        self.maxLines=maxLines                              # 1e10 for everything

        self.persistenceDir=persistenceDir
        

        
        #self.depthDistribution=[0]*self.maxDepthDistribution
        #self.unadjustedpValues={}

        # check whether persistenceDir exists
        if not(os.path.exists(self.persistenceDir)):
            raise IOError("persistenceDir must exist: {0}".format(self.persistenceDir))
        
        self._refresh()
    def _refresh(self):
        """ define data structures used during the parse operation, and connectors """
        self.sampleId=None
        self.vcfSqlObj=None
        self.timings={}
        self.baseproperties={}       
    def _toPandas(self, DictOfDicts):
        """ converts the DictOfDicts into a pandas dataframe """
        df=pd.DataFrame.from_dict(DictOfDicts, orient='index')
        return(df.reindex_axis(sorted(df.columns), axis=1))

    def summaryStatistics(self):
        """ outputs mean, sd, number and median of variant sites """
        x=array(self.sigMinorAlleleFrequencies())
        return {'summary':{'median':np.median(x),'mean':np.mean(x),'std':np.std(x), 'nvariants':len(x)}}
    
    def alleleDictionary(self):
        """ returns a dictionary of bases containing minor allele frequencies.
        The key is the position, and the dictionary contains read counts, and maf.
        
        # might be better as a dataframe?
        
        Suitable for export using DictWriter"""    
        maf={}
        for pos in self.baseproperties.keys():
            
            if self.baseproperties[pos]['qvalue']==0:
                mlp= 250	# if q value recorded as 0, register minus log p as 250
            else:
                mlp= -math.log(self.baseproperties[pos]['qvalue'],10)
                maf[pos]={'base':pos,\
                          'ma_reads':self.baseproperties[pos]['freq_2'],\
                    'depth':self.baseproperties[pos]['depth'],\
                    'maf':float(self.baseproperties[pos]['freq_2'])/float(self.baseproperties[pos]['depth']),\
                    'minusLogQvalue':mlp,\
                    'base_a':self.baseproperties[pos]['base_a'],\
                    'base_c':self.baseproperties[pos]['base_c'],\
                    'base_g':self.baseproperties[pos]['base_g'],\
                    'base_t':self.baseproperties[pos]['base_t']}
        return maf   
    def adjustedpValues(self,x):
        """
        Assumes a list or numpy array x which contains p-values for multiple tests
        Translated from p.adjust function from R.
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
    def _parse(self, vcffile):
        """ parses a vcffile into an in-memory dictionary, self.baseproperties """
        # read and parse the vcf file
        self.baseproperties={}
        with gzip.open(vcffile, "r") as f:
            
            # iterate over the vcf file
            returnedLines=0
            for line in f:
                if line[0] == "#":
                    continue  # it is a comment; go to next line;
                if "INDEL" in line:
                    continue  #this is not needed because ours don't contain INDELs anymore; go to next line;
                
                # parse the line.
                chrom, pos, varID, ref, alts, score, filterx, infos, fields, sampleInfo = line.strip().split()
                pos = int(pos)
                alts = alts.split(",")
                infos = dict(item.split("=") for item in infos.split(";"))
                baseCounts4=map(int, infos['BaseCounts4'].split(","))   #get frequencies of high quality bases
                baseFreqs=map(int, infos['BaseCounts4'].split(","))     #get frequencies of high quality bases
                baseFreqs.sort(reverse=True)                            #get frequencies of high quality bases, sorted
                depth = sum(baseCounts4)
     
                # compute probability from exact binomial test
                if (baseFreqs[0]<depth and depth>0):        # the majority base is not the only base AND depth is more than 0;
                    pvalue=stats.binom_test(x=baseFreqs[1],n=depth,p=self.expectedErrorRate)   # do the test if any variation
                elif baseFreqs[0]==depth:
                    pvalue=1        # there is only one base
                elif depth==0:
                    pvalue=None     # can't tell, no data
                else:
                    raise Error("Logical error: should never reach this point {0} {1}".format(baseFreqs[0], depth))
                
                if pvalue==0:
                    mlp= 250        # code minus log p as 250
                elif pvalue is not None:
                    mlp= -math.log(pvalue,10)
                elif pvalue is None:
                    mlp=None
                    
                # store output in a dictionary   
                if depth>0:
                    maf=float(baseFreqs[1])/float(depth)
                else:
                    maf=None
                self.baseproperties[pos]={'pos':pos, 'ref':ref, 'depth':depth,\
                        'base_a':baseCounts4[0], 'base_c':baseCounts4[1], 'base_g':baseCounts4[2], 'base_t':baseCounts4[3], \
                        'maf':maf,'pvalue':pvalue, 'mlp':mlp}
                            
                returnedLines=returnedLines+1
                if (returnedLines>=self.maxLines):
                    break       # debug setting; we have returned however many lines we need to do our testing;
                if returnedLines % 100000 ==0:
                    print(returnedLines)
                    
            ## apply fdr     
            positions=self.baseproperties.keys()        # which positions we are analysing
            pvalues=[]                                  # extract the p values into a vector
            for position in positions:                  # for all the positions analysed
                pvalue=self.baseproperties[position]['pvalue']
                if not pvalue is None:
                    pvalues.append(self.baseproperties[position]['pvalue'])    # add the unadjusted value to a list
                    
            adjustedpvalues=self.adjustedpValues(pvalues)   # and apply fdr
            
            # write back qvalues into dictionary
            n=-1
            for position in positions:                  # for all the positions analysed
                n+=1
                if not self.baseproperties[position]['pvalue'] is None:
                    qvalue=adjustedpvalues[n]
                    self.baseproperties[position]['qvalue']=qvalue

                if qvalue==0:
                    mlq= 250        # code minus log p as 250
                elif qvalue is not None:
                    mlq= -math.log(qvalue,10)
                elif qvalue is None:
                    mlq=None
                self.baseproperties[position]['mlq']=mlq
                    
    def store(self, vcffile, sampleId, overwrite=False):
        """ processes the base calls from a VCF file.
        Assumes that a 'baseCounts4' element is present, as computed by the U Oxford MMM pipeline,
        and that contains >Q30 bases.
        
        From this, it computes and persists a data table containing allele frequencies and a significance test for each.
        against the null hypothesis that the minor allele frequency in compatible with the sequencing error rate
        being that expected based on Q30.
        
        vcffile  : the .gz vcf file
        sampleId : identifies the sample, normally the guid
        overwrite: if true, deletes any existing database and rebuilds. 
        """
        # test input
        self.timings['start']=time.strftime("%c")       # monitor time taken
        
        if not os.path.exists(vcffile):
            raise IOError("Input vcf file {0} does not exist.".format(vcffile))

        # test whether a database file containing parsed output already exists
        self.sampleId=sampleId
        self.sampleDbFile=None
        dbfile=os.path.join(self.persistenceDir,"{0}.db".format(self.sampleId))
        if os.path.exists(dbfile):
            if overwrite==False:
                self.dbfile=dbfile
                ## TODO: make a connection
                return
            else:
                os.unlink(dbfile)
                
        # parse the vcf
        self._parse(vcffile=vcffile)
        self.sampleId=sampleId
        self.vcfSQL=vcfSQL(sampleId=sampleId, persistenceDir=self.persistenceDir)
        df=self._toPandas(self.baseproperties)
        self.vcfSQL.populate(df)
      
        self.timings['finished']=time.strftime("%c")       # monitor time taken



################### unit tests ###################
class test_vcfparser_init0a(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.*')):
            os.unlink(filename)
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        
class test_vcfparser_init0b(unittest.TestCase):
    def runTest(self):
        targetdir="nonExistentDirectory"
        with self.assertRaises(IOError):
            vp=vcfStore(persistenceDir=targetdir, maxLines=1000)

class test_vcfparser_parse1(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp') 
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp._parse(vcffile=inputfile)
        df=vp._toPandas(vp.baseproperties)
        self.assertTrue(len(df.index)==1000)
        
class test_vcfparser_store1(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.*')):
            os.unlink(filename)
            
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid1')
        
class test_vcfparser_store2(unittest.TestCase):
    def runTest(self):
        inputfile=os.path.join("..",'testdata','nofile.vcf.gz')
        targetdir=os.path.join('..','unitTest_tmp')
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        with self.assertRaises(IOError):
          vp.store(vcffile=inputfile, sampleId='guid1')

class test_vcfparser_store3(unittest.TestCase):
    def runTest(self):
            
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=os.getcwd())
        vp.store(vcffile=inputfile, sampleId='guid2')
        
class test_vcfSQL_init(unittest.TestCase):
    def runTest(self):
        targetDir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetDir, '*.*')):
            os.unlink(filename)
        self.assertFalse(os.path.exists(os.path.join(targetDir, 's1.db')))
        vc=vcfSQL(sampleId='s1', persistenceDir=targetDir)
        self.assertTrue(os.path.exists(os.path.join(targetDir, 's1.db')))
        
class test_vcfSQL_store(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp') 
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp._parse(vcffile=inputfile)
        df=vp._toPandas(vp.baseproperties)
        vc=vcfSQL(sampleId='s1', persistenceDir=targetdir)
        vc.populate(df)
        (nEntries,)=vc.session.query(func.count(vcfBases.pos)).one()
        self.assertEqual(nEntries,1000)


