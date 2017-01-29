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
#import median, array, mean, std, int 
import time
from csv import DictWriter
import math
import cPickle
import unittest
import logging
import pandas as pd
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, Float, DateTime, Boolean, MetaData, select, func, LargeBinary, Index
from sqlalchemy import create_engine
from sqlalchemy import exc  # exceptions


### beginning of vcfSQL definitions

class vcfSQL():
        """ manages a database backend (default: sqlite) for use for VCF storage.  
        You would not normally use this class directly.  Use vcfStore instead, which is a wrapper around vcfSQL.
        """
        def __init__(self, sampleId, persistenceDir ):
            """ creates a connection to an SQLite database called '{sampleId}.db' in persistenceDir.
            If the database does not exist, it is created.          
            Does not trap for errors if the filename, directory etc are invalid
            """
            self.connString='sqlite:///{0}/{1}.db'.format(persistenceDir,sampleId)
            self.persistenceDir=persistenceDir          
            self.engine = create_engine(self.connString)    # connection string for a database
            
            self.md=MetaData()
            self.vcfBases=Table('vcfBases', self.md,
                        Column('pos', Integer, primary_key=True),
                        Column('ref',String(1)),
                        Column('depth',Integer),
                        Column('base_a',Integer),
                        Column('base_c',Integer),
                        Column('base_g',Integer),
                        Column('base_t',Integer),
                        Column('maf',Float, nullable=True),
                        Column('mlp',Float, nullable=True)
                        )

            # a table to hold any lookup values which may be needed, e.g. ranges of bases to exclude if
            # large numbers are involved and they cannot readily be passed in a query, or need to be persisted
            # for later efficient operation.
            self.lookupBases=Table('tmpBases', self.md,
                        Column('pk',Integer, primary_key=True),
                        Column('baseSetId', Integer, index=True),
                        Column('pos', Integer, index=True)
            )
            
            self.vcfBases.create(self.engine, checkfirst=True)
            self.lookupBases.create(self.engine, checkfirst=True)

        def refreshLookupSet(self):
            """ starts with a fresh lookup set """
            self.lookupBases.drop(self.engine)
            self.lookupBases.create(self.engine, checkfirst=True)
            
        def storeLookupSet(self, identifier, posSet):
            """ stores a set of bases into table, lookupBases, identified by an integer identifier. """
            
            content=[]
            for pos in sorted(posSet):
                content.append({'baseSetId':identifier, 'pos':pos})
                
            self.engine.execute(self.lookupBases.insert(), content)
                             
        def populate(self, df):
            """ appends the data in pandas dataframe df into the table vcfBases.
            
            This primary key (pos) is inserted so
            (i) the data should be added in order
            (ii) there must be no duplicates in pos, or an IntegrityError will be raised.
            
            df must exist and be the right format, e.g. as produced by vcfStore.
            required columns are these:
            {'pos','ref','depth','base_a','base_c','base_g','base_t','maf','mlp'}
            Others will not be written.
            
            Insertion is most efficient if the data is not indexed.
            It is recommended to call create_indices only after all data has been loaded,
            perhaps by multiple calls to .populate().
            """
            
            current_columns=set(df.columns)
            required_columns={'pos','ref','depth','base_a','base_c','base_g','base_t','maf','mlp'}
            extra_columns=current_columns - required_columns
            df.drop(extra_columns, axis=1, inplace=True)
            df.to_sql('vcfBases',index=False, if_exists='append', con=self.engine)
            
        def create_indices(self):
                """ creates indices on tables referenced by self.engine.
                Best done after data loading, to keep speeds up """
                i = Index('depth', self.vcfBases.c.depth)
                i.create(self.engine)
                i = Index('mlp', self.vcfBases.c.mlp)
                i.create(self.engine)
                i = Index('maf', self.vcfBases.c.maf)
                i.create(self.engine)
                
class vcfStore():
    """ methods for allowing rapid access to aspects of information within a vcf file.
    
    * extracting base counts from vcf file;
    * testing whether these are compatible with a prior for sequencing error rate;
    * storing them in an indexed, rapid-access manner
    * allowing access to, and computations on, single samples within this store.
    
    """
    def __init__(self, persistenceDir, maxLines=1e10, expectedErrorRate=0.001, chunkSize=100000):
        """ includeAll                : if false, do not report invariant sites in the output
            expectedErrorRate         : the expected error rate from sequencing; used for negative binomial test; 0.01 corresponds to q30
            maxLines                  : lines to read.  Set to 1000 for unittesting; set to 1e10 for everything
            chunkSize                 : the number of rows to persist in memory before writing to disc.
       """
        
        ## parameters
        self.expectedErrorRate=expectedErrorRate            # 0.01 corresponds to q30
        self.maxLines=maxLines                              # 1e10 for everything
        self.chunkSize=chunkSize
        self.persistenceDir=persistenceDir
                
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
    def _parse(self, vcffile):
        """ parses a vcffile.  yields a dictionary, one per base.
        Do not access this directly.  use the store method instead."""
        # read and parse the vcf file
        self.baseproperties={}
        with gzip.open(vcffile, "r") as f:
            
            # iterate over the vcf file
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
                  
                yield({'pos':pos, 'ref':ref, 'depth':depth,\
                        'base_a':baseCounts4[0], 'base_c':baseCounts4[1], 'base_g':baseCounts4[2], 'base_t':baseCounts4[3], \
                        'maf':maf,'pvalue':pvalue, 'mlp':mlp})                                  
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
        
        self.chunksize: number of DNA bases to store in ram before caching to disc.
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
            if overwrite==True:
                os.unlink(dbfile)
            else:
                # exit this routine as the db already exists
                return(0)
                
        # parse the vcf.  Save to database in batches to keep memory requirement low.
        self.sampleId=sampleId
        self.vcfSQL=vcfSQL(sampleId=sampleId, persistenceDir=self.persistenceDir)

        writeDict={}
        writeKey=0
        returnedLines=0
        for item in self._parse(vcffile=vcffile):
                returnedLines=returnedLines+1
                if (returnedLines>self.maxLines):
                    break       # debug setting; we have returned however many lines we need to do our testing;
       
                writeDict[writeKey]=item
                writeKey+=1
                
                if (writeKey % self.chunkSize)==0:
                        #print("Writing chunk .. at {0}".format(returnedLines))
                        df=self._toPandas(writeDict)
                        self.vcfSQL.populate(df)
                        writeKey=0
                        writeDict={}
        
        # write any extra
        if writeKey>0:
                df=self._toPandas(writeDict)
                self.vcfSQL.populate(df)
                
        # place index
        try:
                self.vcfSQL.create_indices()
        except OperationalError:
                # indices already exist
                pass
        
        # finished
        self.timings['finished']=time.strftime("%c")       # monitor time taken

################### unit tests ###################
class test_vcfparser_init0a(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
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
        for (i,item) in enumerate(vp._parse(vcffile=inputfile)):
                pass
                if i==1000:
                        break
    
class test_vcfparser_store1(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
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


class test_vcfSQL_init(unittest.TestCase):
    def runTest(self):
        targetDir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetDir, '*.db')):
            os.unlink(filename)
        self.assertFalse(os.path.exists(os.path.join(targetDir, 's1.db')))
        vc=vcfSQL(sampleId='s1', persistenceDir=targetDir)
        self.assertTrue(os.path.exists(os.path.join(targetDir, 's1.db')))
        
class test_vcfSQL_store3(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp') 
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
      
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)         
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid1')
 
        vc=vcfSQL(sampleId='guid1', persistenceDir=targetdir)
        vc=vcfSQL(sampleId='guid1', persistenceDir=targetdir)
        q=select([func.count(vc.vcfBases.c.pos)])
        (nEntries,)=vc.engine.execute(q).fetchone()
        self.assertEqual(nEntries,1000)
        
class test_vcfSQL_store4(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp') 
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
      
        for filename in glob.glob(os.path.join(targetdir, '*.*')):
            os.unlink(filename)         
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000, chunkSize=300)
        vp.store(vcffile=inputfile, sampleId='guid1')              # test impact of chunksize
 
        vc=vcfSQL(sampleId='guid1', persistenceDir=targetdir)
        q=select([func.count(vc.vcfBases.c.pos)])
        (nEntries,)=vc.engine.execute(q).fetchone()
        self.assertEqual(nEntries,1000)
        
class test_vcfSQL_store5(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp') 
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
      
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)         
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid3')              # test impact of chunksize
 
        vc=vcfSQL(sampleId='guid1', persistenceDir=targetdir)
        vc.create_indices()

class test_vcfparser_store5(unittest.TestCase):
    def runTest(self):
        print("Running load test")
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=os.getcwd())
        vp.store(vcffile=inputfile, sampleId='guid2')

class test_vcfparser_storeLookupSet1(unittest.TestCase):
    def runTest(self):
        """ tests insertion into lookupVcf """
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)
            
        vc=vcfSQL(sampleId='guid1', persistenceDir=targetdir)
        vc.storeLookupSet(identifier=1, posSet=set([1,2,3,4,5,6,7,8,9,10]))
        
        q=select([func.count(vc.lookupBases.c.pos)])
        (nEntries,)=vc.engine.execute(q).fetchone()
        self.assertEqual(nEntries, 10)
        
class test_vcfparser_storeLookupSet2(unittest.TestCase):
    def runTest(self):
        """ tests refreshing the lookupVcf function """
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)
            
        vc=vcfSQL(sampleId='guid1', persistenceDir=targetdir)
        vc.storeLookupSet(identifier=1, posSet=set([1,2,3,4,5,6,7,8,9,10]))
        
        q=select([func.count(vc.lookupBases.c.pos)])
        (nEntries,)=vc.engine.execute(q).fetchone()
        self.assertEqual(nEntries, 10)
        
        vc.refreshLookupSet()
        (nEntries,)=vc.engine.execute(q).fetchone()
        self.assertEqual(nEntries, 0)
        
        vc.storeLookupSet(identifier=1, posSet=set([1,2,3,4,5]))
        q=select([func.count(vc.lookupBases.c.pos)])
        (nEntries,)=vc.engine.execute(q).fetchone()
        self.assertEqual(nEntries, 5)
              