# illustrate depths computed by python script depthDistributions.py and deposited in a directory
rm(list=ls())

# load libraries
library(ggplot2)
library(reshape)
library(xtable)
library(plyr)
library(fpc)  # may need to be installed as root


# specify root directory
rootdir<-'/mnt/nfs/home/dwyllie/TBMIX-real/real'
outputdir<-file.path(rootdir,'output')

## just Louise files
inputfile<-'~/dev/FOREST/src/louiseTBlist.txt'
louise.guid<-read.csv(inputfile, header=FALSE)
names(louise.guid)[1]<-'guid'

# scan the platecsv directory into a single data frame
print('Loading depth histograms ...')
platecsvs.template<-file.path(rootdir,'depths','*.csv')
platecsvs<-Sys.glob(platecsvs.template)

plateguid<-substr(platecsvs, start=46, stop=81)
union(plateguid,as.character(louise.guid$guid))

import.list <- lapply(platecsvs, read.csv, header=FALSE)

print('Merging')
n.loaded<-0
for (this.df in import.list) {
  n.loaded<-n.loaded+1
  if (n.loaded==1) {
    df<-this.df
  } else {
    df<-rbind(df,this.df)
  }
}
names(df)<-c('sampleId','Depth','nReads','allReads','cumProp')
head(df)

p<-ggplot(df, aes(x=Depth,y=cumProp))
p<-p+geom_point(alpha=0.1)
p<-p+geom_path(aes(group=sampleId))
#p<-p+scale_y_log10()
#p<-p+facet_wrap(~sampleId)
#p<-p+coord_cartesian(xlim=c(0,200))
p

cutoff<-0.3
df$overcutoff<-ifelse(df$cumProp>cutoff,1,0)
to.plot<-subset(df, Depth %in% c(5,10,20,30,40,50,60,90))
to.plot$facetby<-sprintf('Depth cutoff = %3s',to.plot$Depth)
labels<-cast(Depth~., value='overcutoff', data=to.plot, fun.aggregate='mean')
names(labels)[2]<-'propbad'
labels$x<-0.7
labels$y<-200
labels$label<-sprintf("%1.2f\n of samples", labels$propbad)
labels$facetby<-sprintf('Depth cutoff = %3s',labels$Depth)

# plot distribution of bad depths
p<-ggplot(to.plot, aes(x=cumProp))
p<-p+geom_histogram()
p<-p+scale_y_continuous('Number of sequences')
p<-p+scale_x_continuous('Proportion of reads over depth cutoff')
p<-p+facet_wrap(~facetby)
p<-p+geom_vline(xintercept=cutoff, lty=2, colour='red')
p<-p+geom_text(aes(x=x,y=y, label=label), data=labels)
p

ggsave(file.path(outputdir, 'depth_quality.svg'))
saveRDS(df, file.path(outputdir, 'depths.Rds')) # selecting Depth==40 & overcutoff==1 yields bad sequences
saveRDS(labels, file.path(outputdir, 'depth_summary.Rds'))

write.table(df, file=file.path(outputdir, 'depths.csv'), sep='\t', row.names=FALSE, quote=FALSE)

