#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=fastqc_chickadees
#SBATCH --mem=40GB
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=1-15

module load fastqc
module load java
module load trimmomatic

if [ $# -lt 1 ]
    then
        echo "
        Runs pre- and post-trim fastqc for samples in samplesheet. 
        Trimming is done with TrimmomaticSE, using the TruSeq2-SE.fa adapter file.
        Trimming parameters are hard coded.
        Script is designed as an array job for a SLURM cluster.

        [-i] Input sample sheet: one sample per line, with full file path and file-type suffix.
        [-p] Path to fastq directory. Must end with a slash.
        [-t] Number of threads.
        "
    else
        while getopts i:p:t: option
        do
        case "${option}"
        in
        i) samplesheet=${OPTARG};;
        p) fastq_dir=${OPTARG};;
        t) threads=${OPTARG};;
        esac
        done

        #samplesheet="/scratch/Users/pein7187/chickadee_fastq/chickadee_fastq_list_rows.txt"
        #fastq_dir=/scratch/Users/pein7187/chickadee_fastq/

        #code for SLURM array job
        fastq=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
        sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}' | sed -e 's/.fastq.gz//g' -e 's/^.\{45\}//'`
        
        #Pre-trim fastqc
        if [ ! -d ${fastq_dir}pre-trim_fastqc ];then
            mkdir ${fastq_dir}pre-trim_fastqc
        fi

        fastqc -t $threads $fastq --outdir ${fastq_dir}pre-trim_fastqc


        #Trimming
        if [ ! -d ${fastq_dir}trimmed_reads ];then
            mkdir ${fastq_dir}trimmed_reads
        fi

        java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar SE -threads $threads -phred33 \
            $fastq ${fastq_dir}trimmed_reads/${sample}.trimmed.fastq.gz \
            ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/TruSeq2-SE.fa:2:30:10 \
            HEADCROP:15 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36


        #Post-trim fastqc
        if [ ! -d ${fastq_dir}post-trim_fastqc ];then
            mkdir ${fastq_dir}post-trim_fastqc
        fi
        
        fastqc -t $threads $fastq --outdir ${fastq_dir}post-trim_fastqc

fi
