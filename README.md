# Automated-BulkRNAseq
This repository provides scripts to automate bulk RNA-seq analysis.

 ### Usage ###
 1.Customize Samples.txt: Create your own Samples.txt file based on your input fastq files and their corresponding labels for differential gene expression analysis. Each line in Samples.txt should contain the name of a sample and its label separated by a tab, as shown below:  
 ERR188044       YRI  
 ERR188104       YRI  
 ERR188234       YRI  
 ERR188245       GBR  
 ERR188257       GBR  
 ERR188383       GBR  

 2.Run BulkRNAseq.sh: Execute the BulkRNAseq.sh script.  
 3.Enter Samples.txt as the input: Provide Samples.txt as the input when prompted.  

 ### Customization ### 
 You need to set up your own Samples.txt file according to your specific input fastq files and sample labels for differential gene expression analysis.  

 ### Prerequisites ###  

 For the shell script:  
 -Fastqc  
 -Trimmomatic  
 -Hisat2  
 -Samtools  
 -htseq-count  
 -hg38_tran reference file  
 -gencode.v38.annotation.gtf reference hg38  

 For the R script:  
 -DESeq2  
 -ggplot2  
 -pheatmap  
 -gplots  
 -reshape2  
 -plyr  

 ### Keynotes ###  
 1.The shell script is designed to analyze both single-end and paired-end fastq files, so don't forget to comment out redundant lines based on your samples.  
 2.Specify all directories based on your file system path.  
 3.The R script needs to be in the same directory where the .count files exist.  

 ### Contact ### 
 For any questions, please contact mahdi.anvari7@ut.ac.ir.  

