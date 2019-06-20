#!/bin/bash

module load STAR

samplesheet="/scratch/Users/pein7187/chickadee_fastq/chickadee_fastq_list_rows.txt"
sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}' | sed -e 's/.fastq.gz//g' -e 's/^.\{45\}//'`

cpus=$SLURM_CPUS_PER_TASK
genomeDir=/scratch/Users/pein7187/STAR_2pass
ref=/scratch/Users/pein7187/taeGut/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa

STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $ref --sjdbFileChrStartEnd /scratch/Users/pein7187/STAR_mapping_1pass/${sample}SJ.out.tab --runThreadN $cpus
