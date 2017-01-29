#!/usr/bin/env python
 
# finds vcf files to process for mixtures;
# generates a series of .toprocess files, which allow multiple threads to be assigned to each .toprocess file.

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

# use a list of samples which are provided in a flat file.
# they are TB samples from either B'ham or Brighton.

# read these from the flat file, and divide the guids into 15 approx equal sized files.
inputfile="/mnt/microbio/HOMES/dwyllie/data/TBMIX2017/input/QCSummaryPHETB_2017-01-29_00-20-41.txt"
df=pd.read_table(inputfile)
print("Read data frame.")
print(df['Sample ID Guuid'])

# divide the guids into 15 output files 
inputdir='/home/local/GEL/dwyllie/data/TBMIX2017'
nfiles=0
outputfiles=dict()
for i in range(15):
    outputfiles[i]=open(os.path.join(inputdir,str(i)+'.toprocess'),'wb')

for (i, guid) in enumerate(df['Sample ID Guuid']):
    print(i, i % 15, guid)
    outputfiles[i % 15 ].write("{0}\n".format(guid))
    
for i in outputfiles.keys():
    outputfiles[i].close()
print("Finished")
