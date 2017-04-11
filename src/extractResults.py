#!/usr/bin/env python
 
# interrogates multiple mixture files.

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

# read the data about the files to be analysed
nGuids=0
sampleIds=list()
for filename in glob.glob("/home/dwyllie/data/TBMIX/db/*.db"):
	guid=os.path.basename(filename)[0:36]
	sampleIds.append(guid)
	nGuids+=1
print("Reading guids, n= {0}".format(nGuids))


print("Loading deep branch information")
tbd_file="/home/dwyllie/dev/TBMIX/refdata/TBDeepBranch.txt"
coll=pd.read_table(tbd_file)
collsites=list(map(int,coll['position']))


# create a summaryStore object;
basedir='/home/dwyllie/data/TBMIX'
targetDir=os.path.join(basedir,'db')
outputDir="/home/dwyllie/data/TBMIX/output"
dbname='sstat'
test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
test_connstring="sqlite:///%s" % test_path
print(test_connstring)

sstat=summaryStore(db=db1, engineName=test_connstring, dbDir=targetDir)


sampleIds2=sampleIds
mmr1=multiMixtureReader(sampleIds=sampleIds2, persistenceDir=targetDir, this_summaryStore=sstat)     #, this_summaryStore=sstat

mmr1.summarise_byFilter(minDepth=50, minP=-1)
print(mmr1.df)
print("Summary complete")
mmr1.df.to_csv(os.path.join(outputDir, 'depthOver50.csv'))


mmr1.summarise_selection(collsites)
print("Summary complete")
print(mmr1.df)
mmr1.df.to_csv(os.path.join(outputDir, 'collsites.csv'))


mmr1.summarise_byFilter(minDepth=50, minP=8)
print(mmr1.df)
print("Summary complete")
mmr1.df.to_csv(os.path.join(outputDir, 'highVariation.csv'))

# compute summaries of all the files using a filter and for a selection;
mmr1.read_selection(collsites)
print(mmr1.df)
maf=mmr1.df[['sampleId','pos','maf','mlp']]
coll_subset=coll[['lineage','position']]
cs=coll_subset.merge(maf, how='left', left_on='position', right_on='pos')
cs.to_csv(os.path.join(outputDir, 'collsites_individually.csv'))

exit()




# cross tabulate
cs_wide1=pd.crosstab(index=cs['sampleId'], columns=cs['lineage'], values=cs['mlp'], aggfunc=sum)
cs_wide1['Sample ID Guuid']=cs_wide1.index
cs_wide2=pd.crosstab(index=cs['sampleId'], columns=cs['lineage'], values=cs['maf'], aggfunc=sum)
cs_wide2['Sample ID Guuid']=cs_wide2.index

dfSubset=mmr1.df[['Location','Acc Number','Collection Date','Samplename','Sample ID Guuid','Platename','Totalreads']]
dfSubset1=dfSubset.merge(cs_wide1, on='Sample ID Guuid', how='left')
dfSubset2=dfSubset.merge(cs_wide2, on='Sample ID Guuid', how='left')


dfSubset2.to_csv(os.path.join(outputDir, 'type2.csv'))
