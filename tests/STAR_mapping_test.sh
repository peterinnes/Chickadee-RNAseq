#!/bin/bash

module load STAR

rundir=/scratch/Users/pein7187/STAR_mapping_test1pass
mkdir $rundir
cd $rundir

STAR --genomeDir /scratch/Users/pein7187/taeGut --readFilesCommand gunzip -c --readFilesIn /scratch/Users/pein7187/chickadee_fastq/BCCH/Chickadee_CU54780_CGATGT_L003_R1_001.fastq.gz --runThreadN 4
