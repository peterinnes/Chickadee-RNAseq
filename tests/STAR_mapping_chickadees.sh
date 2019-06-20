#!/bin/bash

module load STAR



fastq_list=""
for fastq in /scratch/Users/pein7187/chickadee_fastq/BCCH/*; do
    fastq_list+=$fastq','
done
for fastq in /scratch/Users/pein7187/chickadee_fastq/MOCH/*; do
    fastq_list+=$fastq','
done
for fastq in /scratch/Users/pein7187/chickadee_fastq/MOCH/*; do
    fastq_list+=$fastq','
done
fastq_list=$(echo $fastq_list | sed 's/,$//')

#echo $fastq_list
cat /scratch/Users/pein7187/chickadee_fastq/chickadee_fastq_list_rows.txt | while read line; do

sample=$(sed 's/.fastq.gz//g' | sed 's/^.\{45\}//')
runDir=/scratch/Users/pein7187/STAR_mapping_dirs/$sample
mkdir $runDir
cd $runDir

STAR --genomeDir /scratch/Users/pein7187/taeGut --readFilesCommand gunzip -c --readFilesIn $line --runThreadN 15

done
