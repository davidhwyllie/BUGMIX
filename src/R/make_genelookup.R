# load a genbank file.
# generate a data frame consisting of one row for each base, containing the gene at that base.

library(genoPlotR)

# specify root directory
rootdir<-'/mnt/nfs/home/dwyllie/TBMIX-real/real_allbases'
outputdir<-file.path(rootdir,'output')

# load from genbank
url<-'ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.gbk'
refgenome<-read_dna_seg_from_genbank(file=url, tagsToParse='gene')
genome.length<-4411532

# the data frame now contains one row for each gene.
gene.lookup<-data.frame(base=seq(from=1, to=4411532, by=1))
gene.lookup$gene<-'Intergenic'

# mark each gene
refgenome$label<-ifelse(refgenome$name=='NA', refgenome$synonym, refgenome$name)

print("Marking genes")
for (i in 1:nrow(refgenome)) {
  print(refgenome$label[i])
  gene.lookup$gene[refgenome$start[i]:refgenome$end[i]]<-refgenome$label[i]
}

# mark each intergenic region as part of the gene.
# some intergenic regions are v short and cannot be marked individually
print('Marking intergenic regions')
for (i in 1:(nrow(refgenome)-1)) {
  print(refgenome$label[i])
  featstart<-refgenome$end[i]+1
  featend<-refgenome$start[i+1]-1
  gene.lookup$gene[featstart:featend]<-refgenome$label[i]
}
saveRDS(gene.lookup, file=file.path(outputdir,'genelookup.Rds'))

