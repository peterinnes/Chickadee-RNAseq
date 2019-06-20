#!/bin/bash

module load STAR

samplesheet="/scratch/Users/pein7187/chickadee_fastq/chickadee_trimmed_fastq_list.txt"
fastq=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}' | sed -e 's/.trimmed.fastq.gz/_trimmed/g' -e 's/^.\{54\}//'`
cpus=$SLURM_NTASKS
genomeDir=/scratch/Users/pein7187/BCCH_RefGenome

runDir=/scratch/Users/pein7187/STAR_mapping_2passMode/trimmed/BCCH
#mkdir $runDir
cd $runDir

STAR --twopassMode Basic --genomeDir $genomeDir --readFilesCommand gunzip -c --readFilesIn $fastq --outFileNamePrefix /scratch/Users/pein7187/STAR_mapping_2passMode/trimmed/BCCH/${sample}_bcch -runThreadN $cpus
