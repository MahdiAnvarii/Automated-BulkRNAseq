# Automated-BulkRNAseq
 This repository provides scripts to automate bulk RNA-seq analysis.
 
 ### USAGE ###
 -Costumize your own Samples.txt
 -Run BulkRNAseq.R
 -Enter Samples.txt as the input
 
 ### Costumization ###
 You have to set your own Samples.txt file based on your desired input fastq files and the label that
 each sample has for diffrantial gene expression.
 its a normal text file that each line contains the name of a sample and its label separated by tab just like this:
 ERR188044       YRI
 ERR188104       YRI
 ERR188234       YRI
 ERR188245       GBR
 ERR188257       GBR
 ERR188383       GBR

 ### Prerequairments ###
 
 For the shell script:
 -Fastqc
 -Trimmomatic
 -Hisat2
 -Samtools
 -htseq-count
 -hg38_tran reference file
 -gencode.v38.annotation.gtf refernce hg38

 For the R script:
 -DESeq2
 -ggplot2
 -pheatmap
 -gplots
 -reshape2
 -plyr

 ### Keynotes ###
 -The shell script is designed to analysis both single-end and paired-end fastq files, so dont forget to comment out redundent lines based on your samples
 -Specify all directories based on your file system path
 -The R script needs to be in a same directory which the .count files exist 
