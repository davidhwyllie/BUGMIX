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
tbd_file="/home/dwyllie/dev/TBMIX/refdata/Coll2014_LinSpeSNPs_final.csv"
coll=pd.read_table(tbd_file, sep=',', header=0 )

#  make a dictionary of lineage specific positions.
deep_branches = {}
for i in coll.index:
	lineage=str(coll.get_value(i, 'lineage'))
	pos=int(coll.get_value(i,'position'))

	if lineage not in deep_branches.keys():
		deep_branches[lineage]=set()
	deep_branches[lineage].add(pos)


# Coll2014_LinSpeSNPs_final.csv
# create a summaryStore object;
basedir='/home/dwyllie/data/TBMIX'
targetDir=os.path.join(basedir,'db')
outputDir="/home/dwyllie/data/TBMIX/output"

dbname='tbmix'

test_path="<<DEFAULT>>/%s.db" % dbname         
test_connstring="sqlite:///%s" % test_path
test_connstring="postgresql+psycopg2://tbmix:tbmix@localhost:5432/{0}".format(dbname)
print(test_connstring)

sstat=summaryStore(db=db1, engineName=test_connstring, dbDir=targetDir)
sampleIds2= sampleIds
mmr1=multiMixtureReader(sampleIds=sampleIds2, persistenceDir=targetDir, this_summaryStore=sstat)     #, this_summaryStore=sstat

nLineagesTested=0

for lineage in sorted(deep_branches.keys()):
	print(lineage, len(deep_branches[lineage]))
	mmr1.summarise_selection(selection=deep_branches[lineage], this_selectionDescription=lineage)
	nLineagesTested+=1
	if nLineagesTested==1:
		resdf = mmr1.df
	else:
		resdf = resdf.append(mmr1.df, ignore_index=True)
	### BUG: it only opens one database connection
	### and it gives the same answer for everything.

print("Summary complete")

resdf.to_csv(os.path.join(outputDir, 'collsites.csv'))

