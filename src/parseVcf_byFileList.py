#!/usr/bin/python

# tests whether TB sequences are mixed
import argparse
import csv
import os
import sys
import gzip
import glob
import pysam
from scipy import stats
from numpy import median, array, mean, std
import time
from csv import DictWriter
import math
import cPickle
import datetime
from parseVcf import vcfSQL, vcfStore

### objective is to determine the frequency of mixtures across samples
print "Arguments:",sys.argv
basedir='/home/local/GEL/dwyllie/data/TBMIX2017/'
inputfile=os.path.join(basedir,sys.argv[1])
targetdir=os.path.join(basedir,'db')

print(inputfile)
nFiles=0
with open(inputfile,'rb') as f:
    for guid in f.readlines():
        
        guid=guid.strip()
        vcffile="/mnt/microbio/ndm-hicf/ogre/pipeline_output/{0}/MAPPING/2e6b7bc7-f52c-4649-8538-c984ab3894bb_R00000039/STD/basecalls/{0}_v3.vcf.gz".format(guid)
        
        if os.path.exists(vcffile):
            nFiles+=1
            print("Parsing: %s: %s - %s; start %s" % (os.path.basename(inputfile), nFiles,guid,datetime.datetime.now()))					
            vp=vcfStore(persistenceDir=targetdir)
            vp.store(vcffile=vcffile, sampleId=guid)

            print("Complete: %s: %s - %s; start %s" % (os.path.basename(inputfile), nFiles,guid,datetime.datetime.now()))
        else:
            print("Vcf file %s does not exist" % vcffile)
# 					
print('Complete')
exit(0)
