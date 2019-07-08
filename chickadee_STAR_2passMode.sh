#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=STAR2passMode_chickadee
#SBATCH --mem=40GB
#SBATCH --nodes=1
#SBATCH --array=1-3
#SBATCH --output=/scratch/Users/pein7187/Chickadee-RNAseq/slurm_out_err/chickadee_STAR2passMode_%a.out
#SBATCH --error=/scratch/Users/pein7187/Chickadee-RNAseq/slurm_out_err/chickadee_STAR2passMode_%a.err

set -euxo pipefail

if [ $# -lt 1 ]
    then
        echo "
        Aligns trimmed fastq reads to reference genome using STAR in twopassMode.
        Requires previous STAR indexing of reference genome within genomeDir. 
        Output is sorted BAM. Script is designed for running tasks in parallel on a SLURM cluster.

        [-s] Path to sample sheet
        [-r] Path to runDir
        [-g] Path to reference genomeDir
        [-f] Path to trimmed fastqs
        [-o] Path to output dir for sorted bam files
        [-t] Number of threads
        "

    else
        while getopts s:g:r:f:o:t: option
        do
        case "${option}"
        in
        s) samplesheet=${OPTARG};;
        g) genomeDir=${OPTARG};;
        r) runDir=${OPTARG};;
        f) fastqs_path=${OPTARG};;
        o) sortedbamoutdir=${OPTARG};;
        t) threads=${OPTARG};;        
        esac
        done

        sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
        fastq=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}' | sed -e "s@^@$fastqs_path@" -e 's/$/.trimmed.fastq.gz/'`
        #echo $fastq

        if [ -d $runDir ]; then
            cd $runDir
        else 
            mkdir $runDir;
            cd $runDir
        fi
        
        if [ ! -d $sortedbamoutdir ]; then
            mkdir $sortedbamoutdir
        fi

        module load STAR

        STAR --twopassMode Basic \
            --genomeDir $genomeDir \
            --readFilesCommand gunzip -c \
            --readFilesIn $fastq \
            --outFileNamePrefix $sortedbamoutdir${sample}_BCCH_ \
            --outSAMtype BAM SortedByCoordinate \
            --runThreadN 16 #$threads
    fi
