args <- commandArgs(trailingOnly = TRUE)
SamplesFile <- args[1] 
Samples_data <- read.table(SamplesFile, header = FALSE, stringsAsFactors = FALSE)
labels <- factor(c(Samples_data$V2))
#print(labels)

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(gplots)
library(reshape2)
library(plyr)
options(stringAsFactors=F)
files <- list.files("D:/NGS/Samples/P3/Aligned","*.count")
#print(files)

thecounts <- lapply(files, read.delim, header=F)
thecounts <- do.call(cbind,thecounts)
rownames(thecounts) <- thecounts[,1]
thecounts <- thecounts[,-seq(1,ncol(thecounts),2)]
thecounts <- thecounts[-seq(nrow(thecounts)-4,nrow(thecounts)),] 
colnames(thecounts) <- sub("_sorted.count","",files)
#print(head(thecounts))

gr <- labels
colData <- data.frame(group=gr,type="single-end")
cds <- DESeqDataSetFromMatrix(thecounts,colData,design=~group)
cds <- DESeq(cds)
cnt <- log2(1+counts(cds,normalized=T))

dif <- data.frame(results(cds,c("group",levels(gr)[1],levels(gr)[2])))
dif$padj <- p.adjust(dif$pvalue,method="BH")
dif <- dif[order(dif$padj),]
#print(head(dif))

vp <- ggplot(dif,aes(log2FoldChange,-log10(pvalue),color=log2FoldChange))+geom_point()+theme_bw()
ggsave("volcano_plot.pdf", plot = vp, device = "pdf")
