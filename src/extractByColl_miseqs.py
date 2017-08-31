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

# load the list of guids to be analysed.
inputfile = "/home/dwyllie/data/TBMIX/input/guids.txt"
with open(inputfile, 'rt') as f:
	guids= f.readlines()
	guids = [x.strip() for x in guids]
	print("Recovered {0}".format(len(guids)))

# read the data about the files to be analysed
print("Running guid check ..")
sampleIds = list()
missing = 0
found =0
for guid in guids:
	
	filename = "/home/dwyllie/data/TBMIX/db/{0}.db".format(guid)
	if not os.path.exists(filename):
		missing += 1
	else:
		sampleIds.append(guid)
		found += 1
		#if found > 20:
		#	break
		
print("{0} guids are missing; {1} are present.".format(missing, found))



print("Loading deep branch information")
tbd_file="/home/dwyllie/dev/TBMIX/refdata/TBDeepBranch.txt"
coll=pd.read_table(tbd_file)
collsites=list(map(int,coll['position']))


print("Loading deep branch information for deep branches")
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

# create a summaryStore object;

basedir='/home/dwyllie/data/TBMIX'
targetDir=os.path.join(basedir,'db')
outputDir="/home/dwyllie/data/TBMIX/output"
reportDir="/home/dwyllie/data/TBMIX/reports"

# uses postgres; assumes database mixtures exists with
# role tbmix;
# need to
# GRANT ALL PRIVILEGES ON DATABASE mixtures to tbmix;
# GRANT ALL PRIVILEGES ON DATABASE mixtures TO GROUP tbmix_daemon;   # tbmix_daemon is the owner; tbmix is a login role & is a member of tbmix_daemon.

dbname='mixtures'

test_path="<<DEFAULT>>/%s.db" % dbname         
test_connstring="sqlite:///%s" % test_path
test_connstring="postgresql+psycopg2://tbmix:tbmix@localhost:5432/{0}".format(dbname)
print(test_connstring)

sstat=summaryStore(db=db1, engineName=test_connstring, dbDir=targetDir)


nLineagesTested=0
mmr1=multiMixtureReader(sampleIds=sampleIds, persistenceDir=targetDir, this_summaryStore=sstat)     #, this_summaryStore=sstat

for lineage in sorted(deep_branches.keys()):
	print(lineage, len(deep_branches[lineage]))
	mmr1.summarise_selection(selection=deep_branches[lineage], this_selectionDescription=lineage)
	nLineagesTested+=1
	if nLineagesTested==1:
		resdf = mmr1.df
	else:
		resdf = resdf.append(mmr1.df, ignore_index=True)

print("Lineage summary complete")

resdf.to_csv(os.path.join(reportDir, 'collsites.csv'))



exit()

