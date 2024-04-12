#!/bin/bash

trimmomatic(){
    java -jar /home/mahdi/mydir/Packages/Trimmomatic-0.39/trimmomatic-0.39.jar
}

htseq-count(){
    /home/mahdi/anaconda3/lib/python3.11/site-packages/HTSeq-2.0.5-py3.11-linux-x86_64.egg/EGG-INFO/scripts/htseq-count
}

ReadNamesFromFile() {
    local samplesfile="$1"
    while IFS=$'\t' read -r name label || [[ -n $name ]]; do
        names+=("$name")
        labels+=("$label")
    done < "$samplesfile"
}

read -p "Enter your Samples file: " SamplesFile
ReadNamesFromFile "$SamplesFile"

for name in "${names[@]}"; do
    echo "Processing $name"
    CleanedName=$(echo "${name}" | tr -d '[:space:]')
    read kossher

    # Step 1 : QC - Run Fastqc
    # For Single-end read :
    fastqc "/mnt/d/NGS/Samples/P3/Reads/${CleanedName}.fastq" -o "/mnt/d/NGS/Samples/P3/Reads/"
    # For Paired-end read :
    fastqc "/mnt/d/NGS/Samples/P3/Reads/${CleanedName}_1.fastq" -o "/mnt/d/NGS/Samples/P3/Reads/"
    fastqc "/mnt/d/NGS/Samples/P3/Reads/${CleanedName}_2.fastq" -o "/mnt/d/NGS/Samples/P3/Reads/"

    # Step 2 : Trimming - Run Trimmomatic
    # For Single-end read :
    trimmomatic SE -phred33 "/mnt/d/NGS/Samples/P3/Reads/${CleanedName}.fastq" "/mnt/d/NGS/Samples/P3/Reads/${CleanedName}_trimmed.fastq" ILLUMINACLIP:/home/mahdi/mydir/Packages/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8 LEADING:5 TRAILING:5 MINLEN:20
    fastqc "/mnt/d/NGS/Samples/P3/Reads/${CleanedName}_trimmed.fastq" -o "/mnt/d/NGS/Samples/P3/Reads/"
    # For Paired-end read :
    trimmomatic PE -phred33 "/mnt/d/zigene/P2/Reads/${CleanedName}_1.fastq" "/mnt/d/zigene/P2/Reads/${CleanedName}_2.fastq" "/mnt/d/zigene/P2/Reads/${CleanedName}_forward_paired.fastq" "/mnt/d/zigene/P2/Reads/${CleanedName}_forward_unpaired.fastq" "/mnt/d/zigene/P2/Reads/${CleanedName}_reverse_paired.fastq" "/mnt/d/zigene//P2/Reads/${CleanedName}_reverse_unpaired.fastq" ILLUMINACLIP:/home/mahdi/mydir/Packages/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8 LEADING:5 TRAILING:5 MINLEN:20
    fastqc "/mnt/d/zigene/P2/Reads/${CleanedName}_forward_paired.fastq" -o "/mnt/d/zigene/P2/Reads/"
    fastqc "/mnt/d/zigene/P2/Reads/${CleanedName}_reverse_paired.fastq" -o "/mnt/d/zigene/P2/Reads/"

    # Step 3 : Map to reference - Use Hisat2
    # For Single-end read :
    hisat2 -x /mnt/d/NGS/References/hg38_tran/genome_tran -U "/mnt/d/NGS/Samples/P3/Reads/${CleanedName}_trimmed.fastq" -S "/mnt/d/NGS/Samples/P3/Aligned/${CleanedName}.sam" --dta -p 6
    # For Paired-end read :
    hisat2 -x /mnt/d/NGS/References/hg38_tran/genome_tran -1 "/mnt/d/NGS/Samples/P3/Reads/${CleanedName}_forward_paired.fastq" -2 "/mnt/d/NGS/Samples/P3/Reads/${CleanedName}_reverse_paired.fastq" -S "/mnt/d/NGS/Samples/P3/Aligned/${CleanedName}.sam" --dta -p 6

    # Step 4 : Sort Sam file - Use Samtools
    samtools sort "/mnt/d/NGS/Samples/P3/Aligned/${CleanedName}.sam" -o "/mnt/d/NGS/Samples/P3/Aligned/${CleanedName}_sorted.bam"
    samtools view -h "/mnt/d/NGS/Samples/P3/Aligned/${CleanedName}_sorted.bam" > "/mnt/d/NGS/Samples/P3/Aligned/${CleanedName}_sorted.sam"

    # Step 5 : Counting - Use htseq-count
    htseq-count -i gene_name "/mnt/d/NGS/Samples/P3/Aligned/${CleanedName}_sorted.sam" "/mnt/d/NGS/References/gencode.v38.annotation.gtf" > "/mnt/d/NGS/Samples/P3/Aligned/${CleanedName}_sorted.count"

    echo "------------------------------------------------------------------------------------------------------------"
done

Rscript BulkRNAseq.R "$SamplesFile"
