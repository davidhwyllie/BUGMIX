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
from FisherExact import fisher_exact
import numpy  as np
import time
from csv import DictWriter
import math
import unittest
import logging
import warnings
import pandas as pd
import copy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, Float, DateTime, Boolean, MetaData, select, func, LargeBinary, Index
from sqlalchemy import create_engine, select, func
from sqlalchemy import exc  # exceptions
from sqlalchemy.orm import sessionmaker
from scipy.stats import chi2_contingency


### beginning of vcfSQL definitions
## define classes
db1=declarative_base() # classes mapping to persistent database inherit from this
class compare_two_bases():
        """ a class for comparing two base frequencies using a Fisher's exact test """
        def compare(self, dict1, dict2):
                """ compares the frequencies of bases 'A','C','G','T' in each dictionary using a
                Fisher exact test, with code ported from R to Python
                https://github.com/maclandrol/FisherExact/blob/master/README.md
                """
                a1=[dict1['base_a'],dict1['base_c'],dict1['base_t'],dict1['base_g']]
                a2=[dict2['base_a'],dict2['base_c'],dict2['base_t'],dict2['base_g']]
                p=fisher_exact(table=[a1,a2], hybrid=True)
                return(p)
        
class summaryStore():
        """ a class for storing the minor allele frequencies and other statistics from single or groups of samples,
        or genes within them. """        
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
         
            dbDir  = a directory used for creating SQLite database, if specified.
       
            example usage:
            
            dbDir='/home/dwyllie/persistence'
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
            
        def store(self, selectionDescription, mean_maf, mean_mlp, nBases, total_depth, total_nonmajor_depth, mean_depth=None):
                """ stores the selection, provided that selectionDescription does not already exist """
                nBases = int(nBases)
                total_depth = int(total_depth)
                total_nonmajor_depth = total_nonmajor_depth
                this_sstat = summaryStatistics(
                                             selectionDescription=selectionDescription,
                                             mean_maf=mean_maf,
                                             mean_mlp=mean_mlp,
                                             mean_depth=mean_depth,
                                             nBases=nBases,
                                             total_depth=total_depth, 
                                             total_nonmajor_depth=total_nonmajor_depth
                        )
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
                                if not key in ['guid','selectionDescription','mean_depth','nBases','mean_maf','mean_mlp', 'total_depth','total_nonmajor_depth']:
                                        del retVal[key]

                        return(retVal)
        def store_genes(self, df):
                """ stores the dataframe df  """
                print("Storing ..")
                required_fields=set(['gene','sampleId','base_selector','mean_depth','min_depth','max_depth','start','stop','length','mean_mlp','mean_maf', 'total_depth', 'total_nonmajor_depth'])
                if set(df.columns)==required_fields:
                        this_base_selector=min(df.loc[:,'base_selector'].values)
                        this_sampleId=     min(df.loc[:,'sampleId'].values)
                        
                        #print("Storing data for {0} ({1})".format(this_sampleId, this_base_selector))
                        nExisting,=self.session.query(func.count(geneSummaryStatistics.gstatId)).filter(geneSummaryStatistics.guid==this_sampleId).\
                                               filter(geneSummaryStatistics.base_selector==this_base_selector).one()
                        if nExisting==0:
                                for ix in df.index:
                                        print(ix)
                                        if df.loc[ix,'length']>0:
                                                gss=geneSummaryStatistics(
                                                                gene=df.loc[ix,'gene'],
                                                                base_selector=this_base_selector,   
                                                                guid=this_sampleId,    
                                                                mean_depth=int(df.loc[ix,'mean_depth']),
                                                                min_depth=int(df.loc[ix,'min_depth']),
                                                                max_depth=int(df.loc[ix, 'max_depth']),
                                                                start=int(df.loc[ix, 'start']),   
                                                                stop=int(df.loc[ix, 'stop']),   
                                                                length=int(df.loc[ix, 'length']),
                                                                mean_mlp=float(df.loc[ix, 'mean_mlp']),
                                                                mean_maf=df.loc[ix, 'mean_maf'],
                                                                total_depth=int(df.loc[ix, 'total_depth']),
                                                                total_nonmajor_depth=int(df.loc[ix, 'total_nonmajor_depth'])
                                                )
                                                self.session.add(gss)
                                        self.session.commit()

                        else:
                                #print("Records are already present; not stored")
                                pass
                else:        
                        raise TypeError("Data frame must contain the following fields {0} but it looks like this:\n{1}".format(required_fields, df))
            
                return(0)       
        def recover_genes(self, this_base_selector, this_sampleId):
                """ recovers the stored values """
                q=self.session.query(
                        geneSummaryStatistics.gene,\
                        geneSummaryStatistics.base_selector,
                        geneSummaryStatistics.guid,\
                        geneSummaryStatistics.mean_depth,\
                        geneSummaryStatistics.min_depth,\
                        geneSummaryStatistics.max_depth,\
                        geneSummaryStatistics.start,\
                        geneSummaryStatistics.stop,\
                        geneSummaryStatistics.length,\
                        geneSummaryStatistics.mean_mlp,\
                        geneSummaryStatistics.mean_maf,\
                        geneSummaryStatistics.total_depth,\
                        geneSummaryStatistics.total_nonmajor_depth).\
                        filter(geneSummaryStatistics.guid==this_sampleId).\
                        filter(geneSummaryStatistics.base_selector==this_base_selector)
                df=pd.read_sql(sql=q.statement, con=self.engine)
                if len(df.index)==0:
                        return  None
                else:
                        return (df)
        def restart(self):
                """ removes all entries from the database """
                self.session.query(summaryStatistics).delete()
                
                
class summaryStatistics(db1):       # the guids
    """ SQLalchemy class definition for a table including stored summary values """
    __tablename__ = 'mixtureStatistics'
    sstatId=Column(Integer, primary_key=True)
    selectionDescription=Column(String(80), unique=True, index=True)        
    mean_maf=Column(Float, nullable=False)
    mean_mlp=Column(Float, nullable=False)
    mean_depth=Column(Float, nullable=False)
    nBases=Column(Integer, nullable=False)
    total_depth=Column(Integer, nullable=False)    
    total_nonmajor_depth=Column(Integer, nullable=False)   


class geneSummaryStatistics(db1):       # the guids
    """ SQLalchemy class definition for a table including stored summary values """
    __tablename__ = 'perGeneStatistics'
    gstatId=Column(Integer, primary_key=True)
    gene=Column(String(50), index=True)
    base_selector=Column(String(12), index=True)   
    guid=Column(String(36), index=True)    
    mean_depth=Column(Float, nullable=False)
    min_depth=Column(Float, nullable=False)    
    max_depth=Column(Float, nullable=False)
    start=Column(Integer, nullable=False)    
    stop=Column(Integer, nullable=False)    
    length=Column(Integer, nullable=False)
    mean_mlp=Column(Float, nullable=True)
    mean_maf=Column(Float, nullable=True)
    total_depth=Column(Integer, nullable=False)
    total_nonmajor_depth=Column(Integer, nullable=False)
    
class multiMixtureReader():
        """ a class reading the results of mixture computations from multiple TB sequences, as generated by vcfSQL """
        def __init__(self, sampleIds, persistenceDir, this_summaryStore=None, base2geneLookup=os.path.join('..','refdata','refgenome_lookup.txt')):
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
                # instantiate base frequency comparer object
                self.ctb=compare_two_bases()
                
                # load base2gene lookup table, if available
                if base2geneLookup is not None:
                        self.base2gene=pd.read_table(base2geneLookup, engine='c', sep='\t',  compression=None, index_col=0)
                else:
                        self.base2gene=None
                        
                # make a summary store
                if this_summaryStore is None:           # use a database called 'sstat' in the data directory.
                        dbname='sstat'
                        test_path="<<DEFAULT>>/%s.db" % dbname          
                        connstring="sqlite:///%s" % test_path
                        self.summaryStore=summaryStore(db=db1, engineName=connstring, dbDir=persistenceDir)
                else:
                        self.summaryStore=this_summaryStore
                        
                              
                self.mixtureReader={}
                print("Instantiating.  Adding mixtureReaders to the multiMixtureReader")
                for (i,sampleId) in enumerate(sampleIds):
                        if (i % 100) == 0:
                                print("{0} / {1}".format(i,len(sampleIds)))
                        try:

                                single_mmr= mixtureReader(sampleId, persistenceDir, self.summaryStore)
                                self.mixtureReader[sampleId]=copy.copy(single_mmr)
                        except KeyError as e:
                                warnings.warn("no guid {0}".format(sampleId))
                                
                self.df=None

        def report_by_gene(self,guid, report_to_directory=None):
                """ produces a per gene report on mixtures
                
                If report_to_directory is None, stores the result to db (slow, not recommended)"""
                if self.base2gene is None:
                        warnings.warn("No report generated as no base2gene lookup file.")
                        return None     # we have no base2gene lookup table with which to compute per-gene information
                
                if not guid in self.mixtureReader.keys():
                        raise KeyError("{0} has not been loaded.  Instantiate the multiMixtureReader with the necessary guids needed.".format(guid))

                # check whether we have already computed this, either in the database or on disc.
                if report_to_directory is None:
                        # recover from database
                        df=self.summaryStore.recover_genes(this_base_selector='all', this_sampleId=guid)
                        if df is not None:
                                print("stored results found..")              
                                return(df)
                else:
                        outputfile = os.path.join(report_to_directory, 'perGene_{0}.tsv'.format(guid))
                        if os.path.exists(outputfile):
                                df = pd.read_csv(outputfile, sep='\t', index_col=None, header=0)
                                print("stored results found..")
                                return(df)
                        

                # otherwise, we compute the result.
                print('No stored result found.  Reporting by gene for guid {0}'.format(guid))
                self.mixtureReader[guid].read_all()
                l1 = len(self.mixtureReader[guid].df.index)
                l2 = len(self.base2gene.index)
                if not l1 == l2:
                        warnings.warn("Skipping; likely corruption; The base2gene lookup has length {1} but the data set loaded has length {0}".format(l1,l2))
                else:
                        
                        df=pd.concat([self.base2gene.reset_index(drop=True),self.mixtureReader[guid].df], axis=1)              # in R,  this is a cbind operation
                        self.mixtureReader[guid].df=None                                                                       # release memory
         
                        # the process of producing these statistics is now ultrafast 
                        # can readily produce comprehensive per-gene metrics
                        r0= df.groupby(['gene'])['gene'].min().to_frame(name='gene')
                        r0['sampleId']=guid
                        r0['base_selector']='all'
                        
                        r1= df.groupby(['gene'])['depth'].mean().to_frame(name='mean_depth')
                        r2= df.groupby(['gene'])['depth'].min().to_frame(name='min_depth')
                        r3= df.groupby(['gene'])['depth'].max().to_frame(name='max_depth')
                       
                        r4= df.groupby(['gene'])['pos'].min().to_frame(name='start')                
                        r5= df.groupby(['gene'])['pos'].max().to_frame(name='stop')                
                        r6= df.groupby(['gene'])['pos'].count().to_frame(name='length')
                        
                        r7= df.groupby(['gene'])['mlp'].mean().to_frame(name='mean_mlp')
                        
                        # if all mafs are NA, then mean() will fail with a pandas.core.base.DataError
                        try:                
                                r8= df.groupby(['gene'])['maf'].mean().to_frame(name='mean_maf')     
                        except pd.core.base.DataError:
                                r8=r1.copy()
                                r8.columns=['mean_maf']
                                r8['mean_maf']=None

                        # compute total depth
                        r9= df.groupby(['gene'])['depth'].sum().to_frame(name='total_depth')
                        
                        # compute total_nonmajor_depth
                        df['most_common'] = df[['base_a','base_c','base_g', 'base_t']].max(axis=1)
                        df['nonmajor'] = df['depth'] - df['most_common']
                        r10= df.groupby(['gene'])['nonmajor'].sum().to_frame(name='total_nonmajor_depth')
                                                
                        df=pd.concat([r0,r1,r2,r3,r4,r5,r6,r7,r8, r9, r10], axis=1)              # in R,  this is a cbind operation
            
                        print("Computation complete. Storing ..")
                        
                        if report_to_directory is None:
                                self.summaryStore.store_genes(df)
                        else:
                                outputfile = os.path.join(report_to_directory, 'perGene_{0}.tsv'.format(guid))
                                df.to_csv(outputfile, sep='\t', index=False)
                                
                        print("Storing complete.")
                        # *****
        def compare_base_frequencies(self,guid1, guid2, selection):
                """ compares guid1 and guid2's base frequencies  at the positions in the set 'selection'
                using a fisher exact test.
                Returns a vector of log10(p) values """
                if not guid1 in self.mixtureReader.keys():
                        raise KeyError("{0} has not been loaded.  Instantiate the multiMixtureReader with all guids needed.".format(guid1))
                if not guid2 in self.mixtureReader.keys():
                        raise KeyError("{0} has not been loaded.  Instantiate the multiMixtureReader with all guids needed.".format(guid2))
                if not type(selection) in (set, list):
                        raise TypeError("bases must be either a set of a list, but it is a {0}".format(type(bases)))
                
                # read the frequencies into data frames
                self.mixtureReader[guid1].read_selection(selection)
                self.mixtureReader[guid2].read_selection(selection)
                
                pvalues=[]
                mafs=[]
                for i in self.mixtureReader[guid1].df.index:
                        d1=self.mixtureReader[guid1].df.loc[i,].to_dict()
                        d2=self.mixtureReader[guid2].df.loc[i,].to_dict()
                        p=self.ctb.compare(d1,d2)
                        pvalues.append(p)
                        
                        # also store the maximal maf of the pair
                        mafs.append(max(d1['maf'],d2['maf']))
    
                # derive summary statistics
              
                def compute_mlp(x):
                        """ computes -log10(x).  Returns 50 if x<1e-50"""
                        if x<1e-50:
                                return(50)
                        elif x==1:
                                return(0)
                        else:
                                return(-1*math.log10(x))
                        
                retVal={'mafs':mafs,'maf_avg':np.mean(mafs),'maf_median':np.median(mafs)}
                retVal['mlp']=list(map(compute_mlp, pvalues))
                retVal['mlp_sum']=sum(retVal['mlp'])
                retVal['mlp_avg']=np.mean(retVal['mlp'])
                retVal['mlp_median']=np.median(retVal['mlp'])
                
                # which exceed bonferroni adjustment
                p_cutoff=0.01/len(retVal['mlp'])
                retVal['n_over_cutoff']=sum(1 for x in pvalues if x < p_cutoff)
                
                return(retVal)
                
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
                
                # the selection must be integers, not numpy.int64 etc.
                selection=list( map( int, selection ))
                for (i,sampleId) in enumerate(self.mixtureReader.keys()):
                        self.mixtureReader[sampleId].read_selection(selection)
                self._aggregate()
                
        def summarise_selection(self, selection, this_selectionDescription=None):
                """ performs summarise_selection (see mixtureReader) for all sampleIds.
                
                Selection is a set of bases to summarise.
                selectionDescription is an optional, user-friendly description of what the result is about.
                If present, it should be unique to this search, because it will be used to search for cached result(s) to speed computation.
                
                results are put in self.df  """
                
                #print("Summarising")
                for (i,sampleId) in enumerate(self.mixtureReader.keys()):
                        if (i % 100) == 0 & i>0:
                                print("{0} / {1}".format(i,len(self.mixtureReader.keys())))

                        if this_selectionDescription is None:
                                use_selectionDescription = None
                        else:
                           use_selectionDescription = "{0}|{1}".format(sampleId, this_selectionDescription)
                        self.mixtureReader[sampleId].summarise_selection(selection, selectionDescription = use_selectionDescription)

                self._aggregate()
                
        def summarise_byFilter(self, minDepth, minP, this_selectionDescription=None):
                """ performs summarise_byFilter (see mixtureReader) for all sampleIds.               
                minDepth and minP are criteria to select bases for summary.
                selectionDescription is an optional, user-friendly description of what the result is about.
                If present, it should be unique to this search, because it will be used to search for cached result(s) to speed computation.               
                results are put in self.df
                """
                #print("Summarising using a filter")
                for (i,sampleId) in enumerate(self.mixtureReader.keys()):
                        if (i % 100) == 0:
                                print("{0} / {1}".format(i,len(self.mixtureReader.keys())))
                        self.mixtureReader[sampleId].summarise_byFilter(minDepth, minP, this_selectionDescription=this_selectionDescription)
                #print("Aggregating ..")
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

                self.df=None
                               
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
                        
                # positions as integers, or they won't json convert;
                selection=list(map(int, selection))
                
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
                     # then we must compute a selectionDescription from the criteria passed to us.  uses a hash to generate a unique value given inputs
                     selDes=self._generateSelectionDescription(sequenceGuid=self.sampleId, selection=selection, minP=None, minDepth=None)
                else:
                     selDes=selectionDescription                # we use what we're told.  it's the user's responsibility to make sure this is unique.
                     
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
                                nBases=self.df.loc[0,'nBases'],\
                                total_depth=int(self.df.loc[0,'total_depth']),
                                total_nonmajor_depth=int(self.df.loc[0,'total_nonmajor_depth']))
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
                                nBases=self.df.loc[0,'nBases'],\
                                total_depth=self.df.loc[0,'total_depth'],
                                total_nonmajor_depth=self.df.loc[0,'total_nonmajor_depth'])
                        
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
                               func.count(self.vcfBases.c.depth).label('nBases'),\
                               func.sum(self.vcfBases.c.depth).label('total_depth'),\
                               func.sum(
                                 self.vcfBases.c.depth - func.max(
                                        self.vcfBases.c.base_a,
                                        self.vcfBases.c.base_c,
                                        self.vcfBases.c.base_t,
                                        self.vcfBases.c.base_g) 
                                ).label('total_nonmajor_depth')\
                             ]).\
                                where(self.vcfBases.c.pos.in_(selection))
                
                # sum(depth-(maf*depth)), sum(depth)
                #print(self.engine)
                #print(sqlCmd)
                self.df=pd.read_sql_query(sql=sqlCmd, con=self.engine)
                #print(self.engine)
                #print(self.df)          # debug
        def summarise_byFilter_fromdb(self, minDepth, minP):
                """ summarises the statistics on the bases which pass selection by minDepth and minP;
                results are stored in self.df"""
                sqlCmd=sqlCmd=select([func.avg(self.vcfBases.c.maf).label('mean_maf'),\
                               func.avg(self.vcfBases.c.depth).label('mean_depth'),\
                               func.avg(self.vcfBases.c.mlp).label('mean_mlp'),\
                               func.count(self.vcfBases.c.depth).label('nBases'),\
                               func.sum(self.vcfBases.c.depth).label('total_depth'),\
                               func.sum(
                                 self.vcfBases.c.depth - func.max(
                                        self.vcfBases.c.base_a,
                                        self.vcfBases.c.base_c,
                                        self.vcfBases.c.base_t,
                                        self.vcfBases.c.base_g) 
                                ).label('total_nonmajor_depth')\
                             ]).\
                                where(self.vcfBases.c.depth>=minDepth).where(self.vcfBases.c.mlp>=minP)
                #print(self.engine)
                #print(sqlCmd)
                self.df=pd.read_sql_query(sql=sqlCmd, con=self.engine)
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
                line=line.decode()
                if line[0] == "#":
                    continue  # it is a comment; go to next line;
                if "INDEL" in line:
                    continue  #this is not needed because ours don't contain INDELs anymore; go to next line;
                
                # parse the line.
                chrom, pos, varID, ref, alts, score, filterx, infos, fields, sampleInfo = line.strip().split()
                pos = int(pos)
                alts = alts.split(",")
                infos = dict(item.split("=") for item in infos.split(";"))
                baseCounts4=list(map(int, infos['BaseCounts4'].split(",")))   #get frequencies of high quality bases
                baseFreqs=list(map(int, infos['BaseCounts4'].split(",")))     #get frequencies of high quality bases
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
class test_summaryStore_init(unittest.TestCase):
    def runTest(self):
        persistenceDir=os.path.join('..','unitTest_tmp')
        dbname='sstat'
        test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
        test_connstring="sqlite:///%s" % test_path
 
        # delete the contents of unitTest_tmp
        for this_file in glob.glob(os.path.join(persistenceDir, '*.*')):
                os.unlink(this_file)
        sstat=summaryStore(db=db1, engineName=test_connstring, dbDir=persistenceDir)

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
        sstat=summaryStore(db=db1, engineName=test_connstring, dbDir=persistenceDir)
        sstat.restart()         # remove any old records
        if not os.path.exists(os.path.join(persistenceDir, 'sstat.db')):
                self.fail("sstat.db was not created")

        sstat.store(selectionDescription='a', mean_mlp=2, mean_maf=0.5, mean_depth=100, nBases=100, total_depth= 1000, total_nonmajor_depth = 200)
        res=sstat.recover(this_selectionDescription='a')
        self.assertEqual(res['mean_maf'],0.5)
        res=sstat.recover(this_selectionDescription='b')
        self.assertTrue(res is None)
        
        sstat.restart()
        res=sstat.recover(this_selectionDescription='a')        
        self.assertTrue(res is None)        

class test_ctb_1(unittest.TestCase):
    def runTest(self):
        """ tests the compare two bases function """
        ctb=compare_two_bases()
        p=ctb.compare({'base_a':10,'base_c':0,'base_t':190,'base_g':0},{'base_a':10,'base_c':0,'base_t':190,'base_g':0})
        self.assertTrue(p==1)        
     
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
        print(mr.df)
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
        print(mr.df)
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

class test_multiMixtureReader_2(unittest.TestCase):
    def runTest(self):
        """ assesses pairwise comparison using Fisher exact test """
        targetdir=os.path.join('..','unitTest_tmp')
        for filename in glob.glob(os.path.join(targetdir, '*.db')):
            os.unlink(filename)    
        inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
        vp1=vcfStore(persistenceDir=targetdir, maxLines=200)
        vp1.store(vcffile=inputfile, sampleId='guid1')
        inputfile=os.path.join("..",'testdata','006578d3-e834-4b4e-9bdd-620cc1f274ed_v3.vcf.gz')
        vp2=vcfStore(persistenceDir=targetdir, maxLines=200)
        vp2.store(vcffile=inputfile, sampleId='guid2')           
        mmr=multiMixtureReader(sampleIds=['guid1','guid2'],persistenceDir=targetdir)
        s=set(range(100))

        res=mmr.compare_base_frequencies(guid1='guid1',guid2='guid2', selection=s)

class test_multiMixtureReader_3(unittest.TestCase):
    def runTest(self):
        """ assesses gene-by-gene summaries  """
        targetdir=os.path.join('..','testdata')
        for filename in glob.glob(os.path.join(targetdir, 'sstat.db')):         # this is the target database used by default
            os.unlink(filename)
        
            
        mmr=multiMixtureReader(sampleIds=['0a2a14b5-3cd6-4533-81df-1cdac9d66373'],persistenceDir=targetdir)
        
        reportdir=os.path.join('..','unitTest_tmp')
        mmr.report_by_gene(guid='0a2a14b5-3cd6-4533-81df-1cdac9d66373', report_to_directory = reportdir)

class test_multiMixtureReader_4(unittest.TestCase):
    def runTest(self):
        """ assesses gene-by-gene summaries  """
        targetdir=os.path.join('..','testdata')
        for filename in glob.glob(os.path.join(targetdir, 'sstat.db')):         # this is the target database used by default
            os.unlink(filename) 
        mmr=multiMixtureReader(sampleIds=['0a1c6e58-6980-4858-8a55-8b59925e0d19','0a2a14b5-3cd6-4533-81df-1cdac9d66373'],persistenceDir=targetdir)
        mmr.report_by_gene(guid='0a1c6e58-6980-4858-8a55-8b59925e0d19')
        mmr.report_by_gene(guid='0a1c6e58-6980-4858-8a55-8b59925e0d19')
        # should only enter once
        mmr.report_by_gene(guid='0a2a14b5-3cd6-4533-81df-1cdac9d66373')        
        # add another
         
class test_multiMixtureReader_6(unittest.TestCase):
    def runTest(self):
        """ assesses gene-by-gene summaries.  This test case failed initially due to low quality data """
        # create a test summaryStore object;
        targetdir=os.path.join('..','testdata') 
        
        # delete the db if it exists
        if os.path.exists(os.path.join(targetdir, 'sstat.db')):
                os.unlink(os.path.join(targetdir, 'sstat.db'))
                
        dbname='sstat'
        test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
        test_connstring="sqlite:///%s" % test_path
        this_sampleId='0a1c6e58-6980-4858-8a55-8b59925e0d19'
        sstat=summaryStore(db=db1, engineName=test_connstring, dbDir=targetdir)
        mmr=multiMixtureReader(sampleIds=[this_sampleId],
                               persistenceDir=targetdir,
                               this_summaryStore=sstat)
        print("Reporting by gene ..")
        mmr.report_by_gene(guid=this_sampleId)

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
        
        
        
        