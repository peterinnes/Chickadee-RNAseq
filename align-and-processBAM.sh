#!/bin/bash

#SBATCH -p short
#SBATCH --job-name=STAR2passMode_chickadee
#SBATCH --mem=40GB
#SBATCH --nodes=1
#SBATCH --array=1-3
#SBATCH --output=/scratch/Users/pein7187/Chickadee-RNAseq/slurm_out_err/chickadee_STAR2passMode_%a.out
#SBATCH --error=/scratch/Users/pein7187/Chickadee-RNAseq/slurm_out_err/chickadee_STAR2passMode_%a.err

set -exou pipefail

if [ $# -lt 1 ]
    then
        echo "
        Aligns trimmed fastq reads to reference genome using STAR in twopassMode.
        Requires previous STAR indexing of reference genome within genomeDir. 
        Output is a fully processed BAM (sorted, rg-added, duplicate-marked, cigar-split). 
        Script is designed for running tasks in parallel on a SLURM cluster.

        [-s] Path to sample sheet, one sample per line, no path prefix or file-type suffix. Sample name only. 
        [-g] Path to ref genome directory (contains STAR index files as well as ref genome fasta)
        [-f] Path to ref genome fasta.
        [-p] Path to trimmed fastqs: single directory containing all trimmed.fastq.gz files.
        [-o] Path to output dir for sorted/processed bam files.
        [-t] Number of threads.
        [-a] Abbreviation of ref genome used, in order to differentiate b/w data aligned to chickadee ref vs zebra finch ref.
        "

    else
        while getopts s:g:f:p:o:t:a: option
        do
        case "${option}"
        in
        s) samplesheet=${OPTARG};;
        g) genomeDir=${OPTARG};;
        f) ref=${OPTARG};;
        p) fastqs_path=${OPTARG};;
        o) sortedbamoutdir=${OPTARG};;
        t) threads=${OPTARG};;        
        a) refAbrv=${OPTARG};;
        esac
        done
        
        #this is code for parallelization of alignment jobs on a cluster
        sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
        fastq=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}' | sed -e "s@^@$fastqs_path@" -e 's/$/.trimmed.fastq.gz/'`

        if [ ! -d $sortedbamoutdir ]; then
            mkdir $sortedbamoutdir
        fi

        module load STAR
        module load gatk
        module load picard
        
        #align fastqs to ref with STAR
        STAR --twopassMode Basic \
            --genomeDir $genomeDir \
            --readFilesCommand gunzip -c \
            --readFilesIn $fastq \
            --outFileNamePrefix $sortedbamoutdir"$sample"_"$refAbrv"_ \
            --outSAMtype BAM SortedByCoordinate \ #output sorted BAM
            --runThreadN $threads

        #Add read groups
        java -jar /opt/picard/2.6.0/picard-2.6.0.jar AddOrReplaceReadGroups \
            I=$sortedbamoutdir"$sample"_"$refAbrv"_Aligned.sortedByCoord.out.bam \
            O=$sortedbamoutdir"$sample"_"$refAbrv".sorted_RGadded.bam \
            RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample SO=coordinate
        rm $sortedbamoutdir"$sample"_"$refAbrv"_Aligned.sortedByCoord.out.bam

        #Mark duplicates
        java -jar /opt/picard/2.6.0/picard-2.6.0.jar MarkDuplicates \
            I=$sortedbamoutdir"$sample"_"$refAbrv".sorted_RGadded.bam \
            O=$sortedbamoutdir"$sample"_"$refAbrv".sorted_RGadded_dupmarked.bam \
            M=$sortedbamoutdir"$sample"_"$refAbrv".dupmarked_metrics.txt
        rm $sortedbamoutdir"$sample"_"$refAbrv".sorted_RGadded.bam

        #Split Cigars
        gatk SplitNCigarReads \
            -R $ref \
            -I $sortedbamoutdir"$sample"_"$refAbrv".sorted_RGadded_dupmarked.bam \
            -O $sortedbamoutdir"$sample"_"$refAbrv".sorted_RGadded_dupmarked_split.bam
        rm $sortedbamoutdir"$sample"_"$refAbrv".sorted_RGadded_dupmarked.bam

        gzip $sortedbamoutdir"$sample"_"$refAbrv".sorted_RGadded_dupmarked_split.bam
        gzip $sortedbamoutdir"$sample"_"$refAbrv".sorted_RGadded_dupmarked_split.bai
fi
