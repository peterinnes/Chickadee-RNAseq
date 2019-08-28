#!/bin/bash

module load STAR

genomeDir=/scratch/Users/pein7187/BCCH_RefGenome_8-28-19
ref=/scratch/Users/pein7187/Chickadee-RNAseq/BCCH_RefGenome_8-28-19/BCCH_pseudochr_rename.fasta
annotations=/scratch/Users/pein7187/Chickadee-RNAseq/BCCH_RefGenome_8-28-19/genes.gff


STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $ref --sjdbGTFfile $annotations --sjdbGTFtagExonParentTranscript Parent --runThreadN 32
