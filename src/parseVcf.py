#!/usr/bin/python

# tests whether VCF sequences are mixed
# includes class managing persistence of VCF extracts

import csv
import os
import sys
import gzip
import shutil
import glob
import pysam
import json
import hashlib
from scipy import stats
import numpy  as np
import time
from csv import DictWriter
import math
import unittest
import logging
import pandas as pd
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, Float, DateTime, Boolean, MetaData, select, func, LargeBinary, Index
from sqlalchemy import create_engine, select, func
from sqlalchemy import exc  # exceptions
from sqlalchemy.orm import sessionmaker


### beginning of vcfSQL definitions
## define classes
db1=declarative_base() # classes mapping to persistent database inherit from this

class summaryStore():
        """ a class for storing the minor allele frequencies and other statistics from single or groups of samples. """        
        def __init__(self, db, engineName, dbDir):
            """ constructor
            
            db = instance of the orm
            engineName = connection string for the database.  examples:
        
                sqlite:///../sqlitedb/v7.db.  The explicit file path has been included.
                sqlite:///<<DEFAULT>>/v7.db.  module will replace <<DEFAULT>> with the relevant location
                sqlite:///:memory: [not recommended here, as heavy re-querying of ew server will occur]
                
            Other databases, including MySql and SQL server, are supported.
            If this module is to be used in a multi-threaded way, with concurrent requests to summaryStore,
            then Sqlite will become unsuitable and MySql or other database will need to be substituted.
         
            persistenceDir  = a directory used for creating SQLite database, if specified.
       
            example usage:
            
            persistenceDir='/home/dwyllie/persistence'
            dbname='sstat'
            test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
            test_connstring="sqlite:///%s" % test_path
            """
    
            engineName=engineName.replace('<<DEFAULT>>', dbDir)
            ## create or open the database using sql alchemy
            
            self.Base = db1                              # instance of the ORM
            self.engine = create_engine(engineName)     # connection string for a database
            self.Base.metadata.create_all(self.engine)  # create the table(s) if they do not exist
            Session = sessionmaker(bind=self.engine)    # class
            self.session=Session()           
        def store(self, selectionDescription, mean_maf, mean_mlp, nBases, mean_depth=None ):
                """ stores the selection, provided that selectionDescription does not already exist """
                this_sstat=summaryStatistics(selectionDescription=selectionDescription, mean_maf=mean_maf, mean_mlp=mean_mlp, mean_depth=mean_depth, nBases=int(nBases))
                try:
                        self.session.add(this_sstat)
                        self.session.commit()
                except exc.IntegrityError:      # selectionDescription exists
                        self.session.rollback()                
                return(0)
        def recover(self, this_selectionDescription):
                """ recovers the stored values """
                res=self.session.query(summaryStatistics).filter(summaryStatistics.selectionDescription==this_selectionDescription).one_or_none()
                if res is None:
                        return(None)
                else:
                        retVal=res.__dict__

                        # remove SQL alchemy specific key value pairs
                        starting_keys=list(retVal.keys())
                        for key in starting_keys:
                                if not key in ['selectionDescription','mean_depth','nBases','mean_maf','mean_mlp']:
                                        del retVal[key]

                        return(retVal)
        def restart(self):
                """ removes all entries from the database """
                self.session.query(summaryStatistics).delete()           
class summaryStatistics(db1):       # the guids
    """ SQLalchemy class definition for a table including stored summary values """
    __tablename__ = 'mixtureStatistics'
    sstatId=Column(Integer, primary_key=True)
    selectionDescription=Column(String(32), index=True, unique=True)        
    mean_maf=Column(Float, nullable=False)
    mean_mlp=Column(Float, nullable=False)
    mean_depth=Column(Float, nullable=False)
    nBases=Column(Integer, nullable=False)
    
class test_summaryStore_init(unittest.TestCase):
    def runTest(self):
        persistenceDir=os.path.join('..','unitTest_tmp')
        dbname='sstat'
        test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
        test_connstring="sqlite:///%s" % test_path
 
        # delete the contents of unitTest_tmp
        for this_file in glob.glob(os.path.join(persistenceDir, '*.*')):
                os.unlink(this_file)
        sstat=summaryStore(db=db1, engineName=test_connstring, persistenceDir=persistenceDir)

        if not os.path.exists(os.path.join(persistenceDir, 'sstat.db')):
                self.fail("sstat.db was not created")
class test_summaryStore_add(unittest.TestCase):
    def runTest(self):
        persistenceDir=os.path.join('..','unitTest_tmp')
        dbname='sstat'
        test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
        test_connstring="sqlite:///%s" % test_path
 
        # delete the contents of unitTest_tmp
        for this_file in glob.glob(os.path.join(persistenceDir, '*.*')):
                os.unlink(this_file)
        sstat=summaryStore(db=db1, engineName=test_connstring, persistenceDir=persistenceDir)
        sstat.restart()         # remove any old records
        if not os.path.exists(os.path.join(persistenceDir, 'sstat.db')):
                self.fail("sstat.db was not created")

        sstat.store(selectionDescription='a', mean_mlp=2, mean_maf=0.5, mean_depth=100, nBases=100)
        res=sstat.recover(this_selectionDescription='a')
        self.assertEqual(res['mean_maf'],0.5)
        res=sstat.recover(this_selectionDescription='b')
        self.assertTrue(res is None)
        
        sstat.restart()
        res=sstat.recover(this_selectionDescription='a')        
        self.assertTrue(res is None)        
class multiMixtureReader():
        """ a class reading the results of mixture computations from multiple TB sequences, as generated by vcfSQL """
        def __init__(self, sampleIds, persistenceDir, this_summaryStore=None):
                """ creates multiple mixtureReader objects, each corresponding to one item in sampleIds.
                    if a summaryStore object is passed, will attempt to use this to obtain stored results, rather than
                    computation from database files.  This speeds operations on hundred or thousands of files.
                                        
                    Examples of use are below:                               
                                # necessary libraries
                                import os
                                import unittest
                                import uuid
                                import inspect
                                import datetime 
                                import hashlib
                                import random
                                import csv
                                import sys
                                import shutil
                                import glob
                                import gzip
                                import pandas as pd
                                
                                from parseVcf import multiMixtureReader, summaryStore, db1
                                
                                # use a list of samples which are provided in a flat file.
                                # they are TB samples from either B'ham or Brighton.
                                
                                # load reference data
                                print("Loading deep branch information")
                                tbd_file="/mnt/microbio/HOMES/dwyllie/dev/TBMIX/refdata/TBDeepBranch.txt"
                                coll=pd.read_table(tbd_file)
                                collsites=coll['position']
                                
                                # read the data about the files to be analysed
                                print("Loading information on TB samples")
                                inputfile="/mnt/microbio/HOMES/dwyllie/data/TBMIX2017/input/QCSummaryPHETB_2017-01-29_00-20-41.txt"
                                df=pd.read_table(inputfile)
                                sampleIds=df['Sample ID Guuid']
                                
                                # create a summaryStore object;
                                basedir='/home/local/GEL/dwyllie/data/TBMIX2017/'
                                targetDir=os.path.join(basedir,'db')
                                outputDir="/mnt/microbio/HOMES/dwyllie/data/TBMIX2017/output"
                                persistenceDir=outputDir
                                dbname='sstat'
                                test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
                                test_connstring="sqlite:///%s" % test_path
                                sstat=summaryStore(db=db1, engineName=test_connstring, persistenceDir=outputDir)
                                
                                # compute summaries of all the files using a filter and for a selection;
                                sampleIds2=sampleIds
                                mmr1=multiMixtureReader(sampleIds=sampleIds2, persistenceDir=targetDir, this_summaryStore=sstat)
                                mmr1.summarise_byFilter(minDepth=50, minP=8)
                                print(mmr1.df)
                                print("Summary complete")
                                
                                mmr1.summarise_selection(collsites)
                                print("Summary complete")
                                print(mmr1.df)
                                exit()
                                
                                print("There are {0} sampleIds.".format(len(sampleIds)))
                                sampleIds1=sampleIds
                                mmr=multiMixtureReader(sampleIds=sampleIds1, persistenceDir=targetDir)
                                
                                mmr.read_selection(collsites)
                                print(mmr.df)
                                maf=mmr.df[['sampleId','pos','maf','mlp']]
                                coll_subset=coll[['lineage','position']]
                                cs=coll_subset.merge(maf, how='left', left_on='position', right_on='pos')
                                print(cs)
                                
                                # cross tabulate
                                cs_wide1=pd.crosstab(index=cs['sampleId'], columns=cs['lineage'], values=cs['mlp'], aggfunc=sum)
                                cs_wide1['Sample ID Guuid']=cs_wide1.index
                                cs_wide2=pd.crosstab(index=cs['sampleId'], columns=cs['lineage'], values=cs['maf'], aggfunc=sum)
                                cs_wide2['Sample ID Guuid']=cs_wide2.index
                                
                                dfSubset=df[['Location','Acc Number','Collection Date','Samplename','Sample ID Guuid','Platename','Totalreads']]
                                dfSubset1=dfSubset.merge(cs_wide1, on='Sample ID Guuid', how='left')
                                dfSubset2=dfSubset.merge(cs_wide2, on='Sample ID Guuid', how='left')
                                
                                dfSubset1.to_csv(os.path.join(outputDir, 'type1.csv'))
                                dfSubset2.to_csv(os.path.join(outputDir, 'type2.csv'))

                """
                
                # make a summary store
                if this_summaryStore is None:
                        dbname='sstat'
                        test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
                        connstring="sqlite:///%s" % test_path
                        self.summaryStore=summaryStore(db=db1, engineName=connstring, dbDir=persistenceDir)
                else:
                        self.summaryStore=this_summaryStore
                        
                        
                self.mixtureReader={}
                print("Instantiating.  Adding mixtureReaders to the multiMixtureReader")
                for (i,sampleId) in enumerate(sampleIds):
                        if (i % 100) == 0:
                                print("{0} / {1}".format(i,len(sampleIds)))
                        self.mixtureReader[sampleId]=mixtureReader(sampleId, persistenceDir, self.summaryStore)
                self.df=None
                
                
        def _aggregate(self):
                """ examines the data frames produced by all mixtureReaders attached, and sets a data frame self.df
                which is the concatenation of all of these, with the sampleId included in each row. """
                nMr=0
                for sampleId in self.mixtureReader.keys():
                        nMr+=1
                        df=self.mixtureReader[sampleId].df
                        df['sampleId']=sampleId
                        #print(df)               # debug
                        if nMr==1:
                                self.df=df
                        else:
                                self.df=self.df.append(df, ignore_index=True)
                                
        def read_selection(self, selection):
                """ performs read_selection (see mixtureReader) for all sampleIds """
                for (i,sampleId) in enumerate(self.mixtureReader.keys()):
                        if i % 100 ==0:
                                print(i)
                        self.mixtureReader[sampleId].read_selection(selection)
                self._aggregate()
                
        def summarise_selection(self, selection, this_selectionDescription=None):
                """ performs summarise_selection (see mixtureReader) for all sampleIds.
                
                Selection is a set of bases to summarise.
                selectionDescription is an optional, user-friendly description of what the result is about.
                If present, it should be unique to this search, because it will be used to search for cached result(s) to speed computation.
                
                results are put in self.df  """
                print("Summarising")
                for (i,sampleId) in enumerate(self.mixtureReader.keys()):
                        if (i % 100) == 0:
                                print("{0} / {1}".format(i,len(self.mixtureReader.keys())))

                        self.mixtureReader[sampleId].summarise_selection(selection, selectionDescription=this_selectionDescription)
                print("Aggregating ..")
                self._aggregate()
                
        def summarise_byFilter(self, minDepth, minP, this_selectionDescription=None):
                """ performs summarise_byFilter (see mixtureReader) for all sampleIds.               
                minDepth and minP are criteria to select bases for summary.
                selectionDescription is an optional, user-friendly description of what the result is about.
                If present, it should be unique to this search, because it will be used to search for cached result(s) to speed computation.               
                results are put in self.df
                """
                print("Summarising using a filter")
                for (i,sampleId) in enumerate(self.mixtureReader.keys()):
                        if (i % 100) == 0:
                                print("{0} / {1}".format(i,len(self.mixtureReader.keys())))
                        self.mixtureReader[sampleId].summarise_byFilter(minDepth, minP, this_selectionDescription=this_selectionDescription)
                print("Aggregating ..")
                self._aggregate()                
class mixtureReader():
        """ a class reading the results of mixture computations from one TB sequence, as generated by vcfSQL """
        def __init__(self, sampleId, persistenceDir, this_summaryStore=None):
                """ creates a connection to an SQLite database called '{sampleId}.db' in persistenceDir.
                    If the database does not exist, raises KeyError.
                    
                    If passed a summaryStore object, it will persist its results in the result store and attempt to used cached results from the
                    result store in preference to ab initio extraction from the database.
                    
                """
                self.sampleId=sampleId
                dbfile='{0}/{1}.db'.format(persistenceDir,sampleId)
                if not os.path.exists(dbfile):
                        raise KeyError("Database corresponding to {0} does not exist at {1}".format(sampleId, dbfile))
                
                self.connString='sqlite:///{0}'.format(dbfile)
                self.persistenceDir=persistenceDir          
                self.engine = create_engine(self.connString)    # connection string for a database
                self.df=None
                self.summaryStore=this_summaryStore
                self.md=MetaData()
                self.vcfBases=Table('vcfBases',
                        self.md,
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
        def read_all(self):
                """ reads all records for this sampleId into memory as a pandas dataframe.  Note that this can run the machine out of memory. """
                sqlCmd=select([self.vcfBases])
                self.df=pd.read_sql_query(sql=sqlCmd, con=self.engine)        
        def read_selection(self, selection):
                """ reads the results at the positions in selection into memory as a pandas dataframe """
                sqlCmd=select([self.vcfBases]).where(self.vcfBases.c.pos.in_(selection))
                self.df=pd.read_sql_query(sql=sqlCmd, con=self.engine)      
        def _generateSelectionDescription(self, sequenceGuid, selection=None, minDepth=None, minP=None):
                """  generates an md5 hash which reflects the values passed to (selection), minDepth, and minP.   """
                
                # at least one of selection, minDepth and minP must be provided.
                if selection is None and minDepth is None and minP is None:
                        raise ValueError("At least one of selection, minDepth and minP must be provided.")
                if selection is None:
                        selection=[]
                if type(selection)==set:
                        selection=list(selection)
                if type(selection)==pd.core.series.Series:
                        selection=list(selection)
                if not type(selection)==list:
                        raise TypeError("Selection must be either a set or a list, but it is a {0}".format(type(selection)))
                else:
                        selection=sorted(selection)
                
                parameters={'sequenceGuid':sequenceGuid, 'selection':selection, 'minDepth':minDepth,'minP':minP}
                to_hash=json.dumps(parameters, sort_keys=True).encode(encoding='utf-8')
                return(hashlib.md5(to_hash).hexdigest())        
        def summarise_selection(self, selection, selectionDescription=None):
                """ summarises statistics on the bases in {selection}; returns a pandas dataframe.
                Will used cached results to do so if it can.
                
                selection is a set of bases to summarise.
                selectionDescription is an optional, user-friendly description of what the result is about.
                If present, it should be unique to this search, because it will be used to search for cached result(s) to speed computation.
                
                """
                if self.summaryStore is None:         # then it has to compute the summary ab initio
                        warnings.warn("No summary store provided")
                        self.summarise_selection_fromdb(selection)              # sets self.df
                        return()

                # we are instructed to use the summaryStore
                if selectionDescription is None:
                     # then we must compute a selectionDescription from the criteria passed to us.
                     selDes=self._generateSelectionDescription(sequenceGuid=self.sampleId, selection=selection, minP=None, minDepth=None)
                else:
                     selDes=selectionDescription
                     
                # now try to recover
                res=self.summaryStore.recover(this_selectionDescription=selDes) 
                if res is None:
                        # nothing is present
                        #print("nothing recovered")
                        self.summarise_selection_fromdb(selection)      # sets self.df
                        # store self.df
                        #print("Storing newly computed material")
                        self.df.loc[0,'selectionDescription']=selDes
                        #print(self.df)

                        self.summaryStore.store(\
                                selectionDescription=self.df.loc[0,'selectionDescription'],\
                                mean_maf=self.df.loc[0,'mean_maf'],\
                                mean_mlp=self.df.loc[0,'mean_mlp'],\
                                mean_depth=self.df.loc[0,'mean_depth'],\
                                nBases=self.df.loc[0,'nBases'])
                else:
                        # coerce res to a pandas data frame
                        #print("recovered {0}; coercing to pandas".format(res))
                        for key in res.keys():
                                res[key]=[res[key]]                             
                        self.df=pd.DataFrame.from_dict(res, orient='columns')
                        #print(self.df)                               
        def summarise_byFilter(self, minDepth, minP, this_selectionDescription=None):
                """ summarises the statistics on the bases which pass selection by minDepth and minP;
                results are stored in self.df
                Will use cached data if it can.
                
                minDepth and minP are criteria to select bases for summary.
                selectionDescription is an optional, user-friendly description of what the result is about.
                If present, it should be unique to this search, because it will be used to search for cached result(s) to speed computation.
                """
                
                if self.summaryStore is None:         # then it has to compute the summary ab initio
                        #print("No summaryStore provided, computing ab initio")
                        self.summarise_byFilter_fromdb(minDepth, minP)
                        return()

                # we are instructed to use the summaryStore
                if this_selectionDescription is None:
                     # then we must compute a selectionDescription from the criteria passed to us.
                     selDes=self._generateSelectionDescription(sequenceGuid=self.sampleId, selection=None, minP=minP, minDepth=minDepth)
                else:
                     selDes=this_selectionDescription
                #print(selDes)     
                # now try to recover
                res=self.summaryStore.recover(this_selectionDescription=selDes)
                #print(res)
                if res is None:
                        # nothing is present
                        #print('-> Computing ')
                        self.summarise_byFilter_fromdb(minDepth, minP)
                        #print(self.df)
                        self.df.loc[0,'selectionDescription']=selDes

                        self.summaryStore.store(\
                                selectionDescription=self.df.loc[0,'selectionDescription'],\
                                mean_maf=self.df.loc[0,'mean_maf'],\
                                mean_mlp=self.df.loc[0,'mean_mlp'],\
                                mean_depth=self.df.loc[0,'mean_depth'],\
                                nBases=self.df.loc[0,'nBases'])
                        
                else:
                        
                        # coerce res to a pandas data frame
                        #print("-> Recovering")
                        for key in res.keys():
                                res[key]=[res[key]]                             
                        self.df=pd.DataFrame.from_dict(res, orient='columns')
                        
        def summarise_selection_fromdb(self, selection):
                """ summarises statistics on the bases in {selection};
                results are stored in self.df.
                Accesses databases for each file to do so. """
                sqlCmd=select([func.avg(self.vcfBases.c.maf).label('mean_maf'),\
                               func.avg(self.vcfBases.c.depth).label('mean_depth'),\
                               func.avg(self.vcfBases.c.mlp).label('mean_mlp'),\
                               func.count(self.vcfBases.c.depth).label('nBases')\
                             ]).\
                                where(self.vcfBases.c.pos.in_(selection))
                #print(self.engine)
                #print(sqlCmd)
                self.df=pd.read_sql_query(sql=sqlCmd, con=self.engine)
                #print(self.df)          # debug
        def summarise_byFilter_fromdb(self, minDepth, minP):
                """ summarises the statistics on the bases which pass selection by minDepth and minP;
                results are stored in self.df"""
                sqlCmd=sqlCmd=select([func.avg(self.vcfBases.c.maf).label('mean_maf'),\
                               func.avg(self.vcfBases.c.depth).label('mean_depth'),\
                               func.avg(self.vcfBases.c.mlp).label('mean_mlp'),\
                               func.count(self.vcfBases.c.depth).label('nBases')\
                             ]).\
                                where(self.vcfBases.c.depth>=minDepth).where(self.vcfBases.c.mlp>=minP)
                #print(self.engine)
                #print(sqlCmd)
                self.df=pd.read_sql_query(sql=sqlCmd, con=self.engine)
                #print(self.df)
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
        #vp=vcfStore(persistenceDir=os.getcwd())
        #vp.store(vcffile=inputfile, sampleId='guid2')

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
              
class test_mixtureReader_1(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)
            
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid1')
        
        mr=mixtureReader(sampleId='guid1',persistenceDir=targetdir)
        mr.read_all()
        self.assertEqual(len(mr.df.index),1000)

class test_mixtureReader_2(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)
            
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid1')
        
        mr=mixtureReader(sampleId='guid1',persistenceDir=targetdir)
        mr.read_selection(selection=[1,2,3,4,5])
        self.assertEqual(len(mr.df.index),5)
        
class test_mixtureReader_3(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)
            
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid1')
        
        mr=mixtureReader(sampleId='guid1',persistenceDir=targetdir)
        mr.summarise_selection(selection=[1,2,3,4,5])
        self.assertEqual(mr.df.iloc[0]['nBases'],5)

class test_mixtureReader_3b(unittest.TestCase):
    def runTest(self):
        #setup
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)
            
        persistenceDir=targetdir
        dbname='sstat'
        test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
        test_connstring="sqlite:///%s" % test_path
        sstat=summaryStore(db=db1, engineName=test_connstring, persistenceDir=persistenceDir)
 
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid1')
        
        #print("Instantiating mixture reader")
        mr=mixtureReader(sampleId='guid1',persistenceDir=targetdir, this_summaryStore=sstat)

        # should be nothing stored
        (n,)=mr.summaryStore.session.query(func.count(summaryStatistics.sstatId)).one()
        self.assertEqual(n,0)
        
        #print("Summarising")
        mr.summarise_selection(selection=[1,2,3,4,5])
        
        #print("Summarising complete")
        self.assertTrue(mr.df is not None)

        # check something has been stored
        (n,)=mr.summaryStore.session.query(func.count(summaryStatistics.sstatId)).one()
        self.assertEqual(n,1)
        #print(mr.df)
        self.assertEqual(mr.df.loc[0,'nBases'],5)
        
        mr.df=None
        mr.summarise_selection(selection=[1,2,3,4,5])
        self.assertEqual(mr.df.loc[0,'nBases'],5)
        
class test_mixtureReader_4(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)
            
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid1')
        
        mr=mixtureReader(sampleId='guid1',persistenceDir=targetdir)
        mr.summarise_byFilter(minDepth=10, minP=2)
        self.assertEqual(mr.df.iloc[0]['nBases'],82)

class test_mixtureReader_4b(unittest.TestCase):
    def runTest(self):
       #setup
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)
            
        persistenceDir=targetdir
        dbname='sstat'
        test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
        test_connstring="sqlite:///%s" % test_path
        sstat=summaryStore(db=db1, engineName=test_connstring, persistenceDir=persistenceDir)
 
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid1')
        
        #print("Instantiating mixture reader")
        mr=mixtureReader(sampleId='guid1',persistenceDir=targetdir, this_summaryStore=sstat)
        
        mr.summarise_byFilter(minDepth=10, minP=2)
        self.assertTrue(mr.df is not None)
        self.assertEqual(mr.df.loc[0,'nBases'],82)
       
        mr.df=None
        mr.summarise_byFilter(minDepth=10, minP=2)
        self.assertTrue(mr.df is not None)
        self.assertEqual(mr.df.loc[0,'nBases'],82)
            
class test_multiMixtureReader_1(unittest.TestCase):
    def runTest(self):
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)
            
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid1')

        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp=vcfStore(persistenceDir=targetdir, maxLines=1000)
        vp.store(vcffile=inputfile, sampleId='guid2')
        
       
        mmr=multiMixtureReader(sampleIds=['guid1','guid2'],persistenceDir=targetdir)
        
        mmr.summarise_byFilter(minDepth=10, minP=2)
        self.assertEqual(len(mmr.df.index),2)
        
        mmr.summarise_selection([1,2,3,4,5])
        self.assertEqual(len(mmr.df.index),2)

        mmr.read_selection([1,2,3,4,5])
        self.assertEqual(len(mmr.df.index),10)

class test_mixtureReader_generateSelectionDescription_1(unittest.TestCase):
    def runTest(self):
        """ tests the generation of a SelectionDescription using a hash of the values passed to it """
        targetdir=os.path.join('..','unitTest_tmp')    
        mr=mixtureReader(sampleId='guid1',persistenceDir=targetdir)
        with self.assertRaises(ValueError):     
                mr._generateSelectionDescription(sequenceGuid='guid1', selection=None, minP=None, minDepth=None)
                
class test_mixtureReader_generateSelectionDescription_2(unittest.TestCase):
    def runTest(self):
        """ tests the generation of a SelectionDescription using a hash of the values passed to it """
        targetdir=os.path.join('..','unitTest_tmp')    
        mr=mixtureReader(sampleId='guid1',persistenceDir=targetdir)
        res1=mr._generateSelectionDescription(sequenceGuid='guid1',selection=set([1,2,3]), minP=None, minDepth=None)          # passed a set
        self.assertEqual(res1,'3f341ecd8ad0e0734ee6a8a16dbb8650')
        res2=mr._generateSelectionDescription(sequenceGuid='guid1',selection=[1,2,3], minP=None, minDepth=None)               # passed a list
        self.assertEqual(res2,'3f341ecd8ad0e0734ee6a8a16dbb8650')
        res3=mr._generateSelectionDescription(sequenceGuid='guid1',selection=[3,2,1], minP=None, minDepth=None)               # passed a list; check order
        self.assertEqual(res3,'3f341ecd8ad0e0734ee6a8a16dbb8650')                
        res4=mr._generateSelectionDescription(sequenceGuid='guid1',selection=[3,2,1], minP=5, minDepth=None)                  # does minP matter
        self.assertEqual(res4,'c9e1d26f388a667d2d38b18538cc3ee8')
        res5=mr._generateSelectionDescription(sequenceGuid='guid1',selection=[3,2,1], minP=None, minDepth=5)                  # does minDepth matter
        self.assertEqual(res5,'72e94a36ef0f47b34c5811cc0271fb63')
        
        
        
        