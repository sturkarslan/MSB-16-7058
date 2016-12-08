library('DESeq2');library("RColorBrewer"); library("gplots");library('ggplot2');library("genefilter");library('pheatmap');library(gridExtra)
# htseq-counts directory
directory <- "/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/results-03/htseq-counts/"
# load annotations from MicrobesOnline
annotations = read.delim("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=882;export=tab", header=T, sep="\t", stringsAsFactors=F)

#
# Collect htseq count files
sample.WT <- grep("WT",list.files(directory),value=TRUE)
sample.744 <- grep("744",list.files(directory),value=TRUE)
# L-WT samples
sample.L1WT <- grep("L1-WT",list.files(directory),value=TRUE)
sample.L2WT <- grep("L2-WT",list.files(directory),value=TRUE)
sample.L3WT <- grep("L3-WT",list.files(directory),value=TRUE)
sample.L4WT <- grep("L4-WT",list.files(directory),value=TRUE)
# LS-WT samples
sample.LS1WT <- grep("LS1-WT",list.files(directory),value=TRUE)
sample.LS2WT <- grep("LS2-WT",list.files(directory),value=TRUE)
sample.LS3WT <- grep("LS3-WT",list.files(directory),value=TRUE)
# L-744samples
sample.L1744 <- grep("L1-744",list.files(directory),value=TRUE)
sample.L2744 <- grep("L2-744",list.files(directory),value=TRUE)
sample.L3744 <- grep("L3-744",list.files(directory),value=TRUE)
sample.L4744 <- grep("L4-744",list.files(directory),value=TRUE)
# LS-744samples
sample.LS1744 <- grep("LS1-744",list.files(directory),value=TRUE)
sample.LS2744 <- grep("LS2-744",list.files(directory),value=TRUE)
sample.LS3744 <- grep("LS3-744",list.files(directory),value=TRUE)
sample.LS4744 <- grep("LS4-744",list.files(directory),value=TRUE)

# Sample annotations for data frame
input.files = c(sample.L1WT, sample.L2WT, sample.L3WT, sample.L4WT,
                  sample.LS1WT, sample.LS2WT, sample.LS3WT,
                  sample.L1744, sample.L2744, sample.L3744, sample.L4744,
                  sample.LS1744, sample.LS2744, sample.LS3744, sample.LS4744)

input.samples = sub("_htseqcounts.txt", "",input.files)
# Sample names
input.names = c(rep("ST1-WT", 2), rep("ST2-WT", 2), rep("ST3-WT", 2), rep("ST4-WT", 1),
                rep("SR1-WT", 2), rep("SR2-WT", 2), rep("SR3-WT", 2),
                rep("ST1-744", 2), rep("ST2-744", 2), rep("ST3-744", 2), rep("ST4-744", 2),
                rep("SR1-744", 2), rep("SR2-744", 2), rep("SR3-744", 2), rep("SR4-744", 2))
# sample conditions
input.conditions = c("ST1", "ST1", "ST2", "ST2", "ST3", "ST3", "ST4",
                     "SR1", "SR1", "SR2", "SR2", "SR3", "SR3",
                     "ST1", "ST1", "ST2", "ST2", "ST3", "ST3", "ST4", "ST4",
                     "SR1", "SR1", "SR2", "SR2", "SR3", "SR3", "SR4", "SR4")
# strains
input.strains = c(rep("WT", 13), rep("744", 16))

# sample transition states
input.state = c(rep("ST", 7), rep("SR", 6), rep("ST", 8), rep("SR", 8))
# time variable for each transition
input.time = c(1,1,3,3,5,5,7,2,2,4,4,6,6,1,1,3,3,5,5,7,7,2,2,4,4,6,6,8,8)
# replicate info
input.replicates = substr(input.samples, 7, 12)
input.replicates = sub("-WT", "WT",input.replicates)
input.replicates = sub("-744", "744",input.replicates)

# create analysis table
input.table = data.frame(sampleName=input.samples, filename=input.files, names= input.names, condition=input.conditions, strains=input.strains, state=input.state, replicates=input.replicates, time=factor(input.time))
# collect read counts from htseq counts
dds.condition <- DESeqDataSetFromHTSeqCount(sampleTable=input.table, directory=directory, design=~1)
# run DESeq
dds <- DESeq(dds.condition)

##### PCA Plot for all samples
data = plotPCA(rld, intgroup=c("condition","strains", "state"), returnData = T)
percentVar <- round(100 * attr(data, "percentVar"))
p <- ggplot(data, aes(PC1, PC2, color=condition, shape=strains))
p + geom_point(size=4) + #geom_text(aes(PC1, PC2, label=strtrim(as.character(condition), 13)), size=2) +
  labs(x=paste0("PC1: ",percentVar[1],"% variance"), y=paste0("PC2: ",percentVar[2],"% variance"), title="PCA Analysis")

##### Trancript abundance plots for all genes
genes = rownames(dds)
filename = "all-genes-expression.pdf"
# initiate pdf
pdf(file=filename, width = 11, height = 8.5)
# create panels to plot multiple plots on the same page
starts = seq(from = 1, to = length(genes) + 15, by = 15)
for(i in starts){
  start1 = i
  stop1 = start1 + 4
  start2 = stop1 + 1
  stop2 = start2 + 4
  start3 = stop2 + 1
  stop3 = start3 + 4
  cat(start1,"-",stop1,", ", start2,"-",stop2,", ",start3,"-",stop3, "\n")
  p1.genes = genes[start1:stop1]
  p2.genes = genes[start2:stop2]
  p3.genes = genes[start3:stop3]
  
  p1.counts = data.frame()
  for(mygene in p1.genes){
    d <- plotCounts(dds, gene=mygene, intgroup=c("names", "strains", "condition", "replicates"), returnData=TRUE)
    mygeneformatted = sub("DVU_", "DVU", mygene)
    a <- annotations[which(annotations$sysName == mygeneformatted),"desc"]
    if(length(a) == 0){
      p1.counts = rbind(p1.counts, cbind(gene=mygene, d, desc=" "))
    } else {
      p1.counts = rbind(p1.counts, cbind(gene=mygene, d, desc=a))
    }
    p1.counts$desc <- as.character(p1.counts$desc)
  }
  p2.counts = data.frame()
  for(mygene in p2.genes){
    d <- plotCounts(dds, gene=mygene, intgroup=c("names", "strains", "condition", "replicates"), returnData=TRUE)
    mygeneformatted = sub("DVU_", "DVU", mygene)
    a <- annotations[which(annotations$sysName == mygeneformatted),"desc"]
    if(length(a) == 0){
      p2.counts = rbind(p2.counts, cbind(gene=mygene, d, desc=" "))
    } else {
      p2.counts = rbind(p2.counts, cbind(gene=mygene, d, desc=a))
    }
    p2.counts$desc <- as.character(p2.counts$desc)
  }
  p3.counts = data.frame()
  for(mygene in p3.genes[!is.na(p3.genes)]){ # make sure the last genes are not NA
    d <- plotCounts(dds, gene=mygene, intgroup=c("names", "strains", "condition", "replicates"), returnData=TRUE)
    mygeneformatted = sub("DVU_", "DVU", mygene)
    a <- annotations[which(annotations$sysName == mygeneformatted),"desc"]
    if(length(a) == 0){
      p3.counts = rbind(p3.counts, cbind(gene=mygene, d, desc=" "))
    } else {
      p3.counts = rbind(p3.counts, cbind(gene=mygene, d, desc=a))
    }
    p3.counts$desc <- as.character(p3.counts$desc)
  }
  
  p1 <- ggplot(p1.counts, aes(x=factor(condition, levels=c("ST1", "SR1","ST2","SR2","ST3","SR3","ST4","SR4")), y=count, color=strains, group=replicates)) +
    geom_point(size=1) +
    scale_y_log10(breaks=c(1,10,100,1000)) +
    geom_line(aes(color=strains), size=1) +
    theme_gray() +
    labs(x="", y="") + 
    theme(axis.text.x=element_text(angle=90, hjust=0,vjust=0, size=7), strip.text.x = element_text(size=8)) +
    facet_grid(.~ gene) + facet_wrap( gene ~ desc, ncol = 5)
  
  p2 <- ggplot(p2.counts, aes(x=factor(condition, levels=c("ST1", "SR1","ST2","SR2","ST3","SR3","ST4","SR4")), y=count, color=strains, group=replicates)) +
    geom_point(size=1) +
    scale_y_log10(breaks=c(1, 10,100,1000)) +
    geom_line(aes(color=strains), size=1) +
    theme_gray() +
    labs(x="", y="") + 
    theme(axis.text.x=element_text(angle=90, hjust=0,vjust=0, size=7), strip.text.x = element_text(size=8)) +
    facet_grid(.~ gene) + facet_wrap( gene ~ desc, ncol = 5)
  
  p3 <- ggplot(p3.counts, aes(x=factor(condition, levels=c("ST1", "SR1","ST2","SR2","ST3","SR3","ST4","SR4")), y=count, color=strains, group=replicates)) +
    geom_point(size=1) +
    scale_y_log10(breaks=c(1,10,100,1000)) +
    geom_line(aes(color=strains), size=1) +
    theme_gray() +
    labs(x="", y="") + 
    theme(axis.text.x=element_text(angle=90, hjust=0,vjust=0, size=7), strip.text.x = element_text(size=8)) +
    facet_grid(.~ gene) + facet_wrap( gene ~ desc, ncol = 5) 
  
  #plot_grid(p1, p2, p3, nrow = 3)
  grid.arrange(p1, p2, p3, nrow=3)
  
}
dev.off()

## Data transformations
# rlog transformation
rld <- rlog(dds)
# variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds)

# Quality checks
# select highly variant top 20 genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
# no transformation
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("condition","strains")])
# normalized log2 transform
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
# rlog transform
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
# variance stabilizing transformation
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)

# Sample to sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$strains, sep="-")
colnames(sampleDistMatrix) <- paste(rld$condition, rld$strains, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, show_colnames=TRUE)