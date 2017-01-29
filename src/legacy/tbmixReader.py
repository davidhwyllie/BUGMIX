#!/usr/bin/env python
import glob
import os
import collections
import csv

# class to read mapping summary data produced by TBMIX.
class mixtureReader():
    """  class to read mapping summary data produced by TBMIX.
    
    # Example usage:
        mr=mixtureReader()
        guid='fff9ea54-56a8-42f2-be10-6e8ad471480b'
        mr.readSummary(guid)
        mr.readDepthDistribution(guid)
        mr.readVariantSites(guid)
        print mr.annotations
            
    No unit tests as yet.
    """
    
    def __init__(self, basedir='/home/local/GEL/dwyllie/data/TBMIX-real/201603/mixedbases'):
        self._refresh()
        self.basedir=basedir    # basedir is the directory in which the mixture files are located, e.g. ~/data/TBMIX-real/201603/tbquant
    def _refresh(self):
        """ sets various properties to the 'factory' (empty) setting """
        self.annotations={}
        self.guid=None
        self.minusLogQvalueFilter=minusLogQvalueFilter
        self.depthFilter=depthFilter
    
    def readSummary(self, guid):
        """ reads a summary file """
        
        if not self.guid==guid:
            self._refresh()            # new guid; remove old data.
            self.guid=guid
            
        searchpath=os.path.join(self.basedir, guid, 'summary.csv')
        inputfiles=glob.glob(searchpath) #inputfile should be here    
        if len(inputfiles)==1:
            # process and return
            with open(inputfiles[0],'rb') as f:
                lines=f.readlines()
                lines=[x.strip('\n') for x in lines]
                if len(lines)==2:
                    tags=lines[0].split(',')
                    linebits=lines[1].split(',')
                    for i in range(len(tags)):
                        self.annotations['mixtures:Summary_'+str(tags[i].strip('\r').strip('\n'))]=linebits[i].strip('\r').strip('\n')                  
                    self.annotations['mixtures:summary_mean']='succeeded'

        else:
            # return no data
            self.annotations['mixtures:readSummary':'failed']
            return(1)
        return(0)
    
    def readDepthDistribution(self, guid):
        """ reads a histogram file """
        
        if not self.guid==guid:
            self._refresh()            # new guid; remove old data.
            self.guid=guid
            
           
        searchpath=os.path.join(self.basedir, guid, 'depthDistribution.csv')
        inputfiles=glob.glob(searchpath) #inputfile should be here    
        if len(inputfiles)==1:
            # process and return
            with open(inputfiles[0],'rb') as f:
                lines=f.readlines()
                lines=[x.strip('\n') for x in lines]
                depths=[]
                totdepth=0
                for i in range(2,len(lines)):   # ignore zero depths
                    line=lines[i].split(',')
                    linebits=lines[i].split(',')
                    depths.append(int(linebits[1]))
                peak=(max(set(depths)))
                for i in range(len(depths)):
                    if depths[i]==peak:
                        self.annotations['mixtures:modalDepth']=i                  
                self.annotations['mixtures:readDepthDistribution']='succeeded'
        else:
            # return no data
            self.annotations['mixtures:readDepthDistribution':'failed']
            return(1)
        return(0)
    
    def readVariantSites(self, guid, minusLogQvalueFilter=10, depthFilter=50): 
        """ reads a variant sites output file .
        Reports mean maf and number of reads passing and not passing the minusLogQ and depth filters set. """
        
        if not self.guid==guid:
            self._refresh()            # new guid; remove old data.
            self.guid=guid
            self.minusLogQvalueFilter=minusLogQvalueFilter
            self.depthFilter=depthFilter
            
        searchpath=os.path.join(self.basedir, guid, 'variantSites.csv')
        inputfiles=glob.glob(searchpath) #inputfile should be here    
        if len(inputfiles)==1:
            # process and return
            with open(inputfiles[0],'rb') as f:
                reader=csv.DictReader(f)
                mafs={'Q+D+':[],'Q-D-':[],'Q+D-':[],'Q-D+':[]}
                for row in reader:
                    if float(row['minusLogQvalue'])>=minusLogQvalueFilter:
                        category='Q+'
                    else:
                        category='Q-'
                    if float(row['depth'])>=depthFilter:
                        category=category+'D+'
                    else:
                        category=category+'D-'
                    mafs[category].append(float(row['maf']))
        
                        
            for key in ['Q+D+']:       #mafs.keys():          # only report double positive bases
                self.annotations['mixtures:'+key+'_nVariantBases']=len(mafs[key])
                if len(mafs[key])>0:
                    self.annotations['mixtures:'+key+'_avgMaf']=sum(mafs[key])/len(mafs[key])
                else:
                    self.annotations['mixtures:'+key+'_avgMaf']=0                  
            self.annotations['mixtures:readVariantSites']='succeeded'
        else:
            # return no data
            self.annotations['mixtures:readVariantSites':'failed']
            return(1)
        return(0)
    

            

        

