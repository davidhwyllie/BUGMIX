# analyse minor species in TB samples
rm(list=ls())
library(ggplot2)
library(reshape)

rootdir<-'/mnt/nfs/home/dwyllie/TBMIX-real/real_allbases'
gldir<-'/mnt/nfs/home/dwyllie/TBMIX-real/real_allbases/output'

# load information about the genome
print('Loading mask file')
genomelength<-4411532

maskfile<-'/mnt/nfs/home/hangphan/refs/R00000039/R00000039_repmask.array.orig'
maskStatus<-read.table(maskfile)
names(maskStatus)[1]<-'isMasked'
maskStatus$base<-1:nrow(maskStatus)
isMasked<-maskStatus$base[maskStatus$isMasked==1]

# read the base to gene lookup table created by make_genelookup.R
print('Loading position to gene lookup ..')
outputdir<-file.path(rootdir,'output')
inputpath1<-file.path(rootdir,'mixedbases','*','variantSites.csv')
gl<-readRDS(file.path(gldir, 'genelookup.Rds'))
gl<-merge(gl,maskStatus, by='base')

basecallsummaries<-Sys.glob(inputpath1)
print(paste("There are ",length(basecallsummaries),"files to analyse"))

# load data, model
print("Scanning .. loading base call summaries for modelling")
vartot<-function(x) {list(nVariants=length(x),maf=mean(x))}

# compute middle of each gene
medl<-function(x) { list(gene.centre=median(x), nBases=length(x))}
middle.of.each.gene<-cast(gene~., value='base',data=gl, fun.aggregate=medl)

isMaskedOptions<-list('All'=c(0,1),'Not masked'=c(0))

# compute summary stats per gene
summary.stats<-function(x) { 
  x<-x[!is.na(x)]
  if (length(x)>1) {
    retval<-list(n=length(x),mean=mean(x), sd=sd(x), upperCI=mean(x)+1.96*(sd(x)/sqrt(length(x))), lowerCI=mean(x)-1.96*(sd(x)/sqrt(length(x))), median=median(x))
  } else {
    retval<-list(n=length(x),mean=mean(x), sd=0, upperCI=mean(x), lowerCI=mean(x), median=median(x))
  }
  retval
}

accessory<-read.csv("~/data/readsvstb/website/allcsv/sampleresults.csv")
accessory<-subset(accessory,speciation.status=='TB')
sampleIds<-unique(as.character(accessory$sampleId))

# iterate over depth cutoffs
n.depths<-0
for (this.depth in c(50)) {  # 0,20,30,40,50,75,100,200
# iterate over inclusion or exclusion of masked

    n.files<-0
    counts<-list()
    for (basecallsummary in basecallsummaries) {
      
      # can put in a depth and maskedFilter here
      test<-scan(basecallsummary, what=character(), quiet=TRUE, nlines=2)
      if (length(test)>0) {
        readcounts<-read.csv(basecallsummary)
        
        readcounts<-subset(readcounts, depth>this.depth)
        readcounts<-merge(readcounts, gl, by='base')
      
        if (nrow(readcounts)>0) {
          tag<-unique(as.character(readcounts$sampleId))
          print(paste(this.depth,n.files, tag))
          
          if (tag %in% sampleIds) { 
            # it's TB
              n.files<-n.files+1
              
              count.per.gene<-cast(gene~., value='maf', fun.aggregate='vartot', data=readcounts)
              count.per.gene<-merge(middle.of.each.gene,count.per.gene, by='gene',all.x=TRUE)
              are.zero<-which(is.na(count.per.gene$maf))
              if (length(are.zero)>0) {
                count.per.gene$nVariants[are.zero]<-0
                count.per.gene$maf[are.zero]<-0
              }      
              count.per.gene$sampleId<-tag
              if (n.files==1) {
                counts.per.gene<-count.per.gene
              } else {
                counts.per.gene<-rbind(counts.per.gene,count.per.gene)
              }
            }
          }
        }
      }
    }
    
    #print('Saving counts per gene')
    #saveRDS(counts.per.gene, file=file.path('/home/dwyllie/data/TBMIX-real/real/fits','countspergene.Rds'))
    if (nrow(counts.per.gene)>0) {
      
      n.depths<-n.depths+1
      
      counts.per.gene$value<-counts.per.gene$nVariants/counts.per.gene$nBases
      counts.per.gene$sampleId<-factor(counts.per.gene$sampleId)
      
      counts.per.gene$value<-counts.per.gene$nVariants/counts.per.gene$nBases
      gene.summary<-cast(gene~., data=counts.per.gene, fun.aggregate="summary.stats")
      gene.summary<-merge(gene.summary,middle.of.each.gene,by='gene')
      gene.summary$depthCutoff<-this.depth
      
      if (n.depths==1) {
        gene.summaries<-gene.summary
      } else {
        gene.summaries<-rbind(gene.summaries, gene.summary)
      }
    }
}

## could convert into an sqlite backed datastore

saveRDS(gene.summaries, file=file.path('/home/dwyllie/data/TBMIX-real/real_allbases/output','genesummary.Rds'))


p<-ggplot(gene.summaries, aes(x=gene.centre, y=mean))
p<-p+geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), data=subset(gene.summaries, upperCI>0))
p<-p+geom_point()
p<-p+scale_x_discrete(name='Chromosome position',breaks=NULL)
p<-p+scale_y_continuous(name='Mean (95% CI) bases mixed per gene')
p
ggsave(p, file=file.path('/home/dwyllie/data/TBMIX-real/real_allbases/output','pergeneplot.png'))


p<-ggplot(gene.summaries, aes(x=mean))
p<-p+geom_density(aes())
p<-p+geom_rug()
p<-p+scale_x_continuous(name='Mean bases mixed per gene')
p
ggsave(p, file=file.path('/home/dwyllie/data/TBMIX-real/real_allbases/output','pergeneplot1.png'))

p<-ggplot(gene.summaries, aes(x=mean))
p<-p+geom_density(aes())
p<-p+geom_rug()
p<-p+scale_x_continuous(name='Mean bases mixed per gene', limits=c(0,0.0025))
p<-p+geom_vline(xintercept=0.001, lty=2, colour='red')
p
ggsave(p, file=file.path('/home/dwyllie/data/TBMIX-real/real_allbases/output','pergeneplot2.png'))

View(gene.summaries)
View(subset(gene.summaries, mean>0.001))

gene.summaries<-gene.summaries[order(-gene.summaries$mean),]

