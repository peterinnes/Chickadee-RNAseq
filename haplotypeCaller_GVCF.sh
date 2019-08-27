#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=call_haplotypes
#SBATCH --mem=40GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-15
#SBATCH --output=/scratch/Users/pein7187/Chickadee-RNAseq/new_pipeline/slurm_out_err/joint_call_haplotypes_%a.out
#SBATCH --error=/scratch/Users/pein7187/Chickadee-RNAseq/new_pipeline/slurm_out_err/joint_call_haplotypes_%a.err

module load picard
module load gatk
module load samtools

if [ $# -lt 1 ]
then
    echo "
    Calls variants with GATK HaplotypeCaller in -GVCF mode.
    Input sorted, processed bam files.
    Settings are per GATK recommendations for RNA-seq data.
    Output is one .gvcf table per sample. Genotyping is done in subsequent step.
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
    r) ref=${OPTARG};; 
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
fi
