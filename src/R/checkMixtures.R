# analyse minor species in TB samples
library(ggplot2)

vcfdir<-'/home/dwyllie/TBMIX/real/vcfcache'
outdir<-'/home/dwyllie/TBMIX/real/vcfplots'
vcfsummaries<-Sys.glob(file.path(vcfdir,'*.txt'))
genomelength<-4411532

reorder<-function(x) {x[order(-x)]}
binom.tests<-function(x) { 
  # if x[1]=0 and x[2]>0, p=1
  if (x[1]==0 & x[2]>0) {
    return(1)
  }
  # otherwise to binomial test 
  if (x[2]>0) { # if there are some reads
    bt<-binom.test(x=x[1],n=x[2], p=x[3],alt='greater')
    retval<-bt$p.value } else 
    { # if there are no reads,return NA
      retval<-NA
    }
  return(retval)
}

msel<-function(stratum,x) {
  c(stratum,mean(x),sd(x),length(x))
}

padtogenomelength<-function(x, genomelength) {
  stopifnot(length(x)<=genomelength)
  x<-c(x,rep(0,genomelength-length(x)))
  x
}

replace.all<-TRUE
maxrows=-1  # set to -1 for all
sequencingerror=1e-3
n.files<-0

for (vcfsummary in vcfsummaries) {
  n.files<-n.files+1
  sampleId<-gsub('.txt','',basename(vcfsummary), fixed=TRUE)
  print(sampleId)
  
  # set the rdata file we are supposed to be making
  rdatafile<- file.path(outdir, paste(sampleId,'.Rds',sep=''))
                        
  # if the rdata file exists, skip and move on.
  if (replace.all==TRUE | !file.exists(rdatafile)) {
    readcounts<-read.csv(vcfsummary, nrows=maxrows)
    names(readcounts)<-c('base','ref','nreads','a','c','g','t')
    readcounts$sampleId<-sampleId
    
    # sort each row
    df<-data.frame(t(apply(readcounts[,4:7], 1, reorder)))
    names(df)<-c('c1','c2','c3','c4')
    readcounts<-cbind(readcounts,df) 
  

    m<-data.frame(c2=readcounts$c2, n=readcounts$nreads, sequencingerror=sequencingerror)

    readcounts$p.value.bt<-apply(m, 1, binom.tests)
    readcounts$q.value<-p.adjust(readcounts$p.value.bt,method='fdr')
    readcounts$obs.prop<-readcounts$c2/readcounts$nreads
    readcounts$odds.minor<-log10(readcounts$c2/(readcounts$nreads-readcounts$c2))
    readcounts$isVariant<-factor(ifelse(readcounts$q.value<0.01,'Variant','Invariant'))    
    
    saveRDS(readcounts, file=rdatafile)
  }   
}

