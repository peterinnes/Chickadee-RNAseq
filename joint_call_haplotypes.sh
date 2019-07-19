#!/bin/bash

module load picard
module load gatk
module load samtools

#old variables before I added commandline arguments
#samplesheet=/scratch/Users/pein7187/chickadee_samplesheet.txt
#sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
#bam_path=/scratch/Users/pein7187/sorted_RGadded_dupmarked_split_bams/zebra_finch/
#bam_suffix=_trimmed_tG.sorted_RGadded_dupmarked_split.bam
#ref=/scratch/Users/pein7187/taeGut/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa
#outdir=/scratch/Users/pein7187/vcf/zebra_finch/joint_calling/

if [ $# -lt 1 ]
then
    echo "
    Calls variants with GATK HaplotypeCaller in -GVCF mode Input sorted, processed bam files.
    Settings are per GATK recommendations for RNA-seq data.
    Output is a VCF table. Genotyping is done in subsequent step.
    Script designed to run on a SLURM cluster.
    
    [-s] Path to samplesheet, one sample per line, no path prefix or file-type suffix. Sample name only.
    [-b] Path to directory containing sorted/processed bams
    [-r] Path to reference fasta
    [-a] abbreviation of reference genome used, to keep track of files aligned to different refs.
    [-o] Path to vcf output directory
    "

else
    while getopts s:b:r:a:o: option
    do
    case "${option}"
    in
    s) samplesheet=${OPTARG};;
    b) bams_path=${OPTARG};;
    r) genomeDir=${OPTARG};; 
    a) refAbrv=${OPTARG};;
    o) vcf_outdir=${OPTARG};;
    esac
    done

#SLURM job array syntax
sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`

gatk HaplotypeCaller \
-R $ref \
-I $bams_path"$sample"_"$refAbrv".sorted_RGadded_dupmarked_split.bam \
-O $vcf_outdir"$sample".gvcf \
-ERC GVCF \
--dont-use-soft-clipped-bases \
--standard-min-confidence-threshold-for-calling 20
