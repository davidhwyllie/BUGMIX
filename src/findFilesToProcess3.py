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

# read these from the flat file, and divide the guids into approx equal sized files.
to_process=[]
globpattern="/mnt/microbio/ndm-hicf/ogre/pipeline_output/*/MAPPING/2e6b7bc7-f52c-4649-8538-c984ab3894bb_R00000039/STD/basecalls/*v3.vcf.gz"
for inputfile in glob.glob(globpattern):
    guid=os.path.basename(inputfile)[0:36]
    # test whether the file has already been parsed
    if os.path.exists(os.path.join('/home/local/GEL/dwyllie/data/TBMIX2017/db','{0}.db'.format(guid))):
        print('exists {0}'.format(guid))
    else:
        print('queuing for processing: {0}'.format(guid))
        to_process.append(guid)

# divide the guids into 10 output files 
inputdir='/home/local/GEL/dwyllie/data/TBMIX2017'
nfiles=0
outputfiles=dict()
for i in range(10):
    outputfiles[i]=open(os.path.join(inputdir,str(i)+'.toprocess'),'wb')

for (i, guid) in enumerate(to_process):
    print(i, i % 10, guid)
    outputfiles[i % 10 ].write("{0}\n".format(guid))
    
for i in outputfiles.keys():
    outputfiles[i].close()
print("Finished")
