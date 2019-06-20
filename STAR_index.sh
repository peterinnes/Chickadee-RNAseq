#!/bin/bash

module load STAR

genomeDir=/scratch/Users/pein7187/BCCH_RefGenome
ref=/scratch/Users/pein7187/BCCH_RefGenome/BCCH_pseudochr_rename.fasta
annotations=/scratch/Users/pein7187/BCCH_RefGenome/snap2.all.gtf

STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $ref --sjdbGTFfile $annotations --runThreadN 32
