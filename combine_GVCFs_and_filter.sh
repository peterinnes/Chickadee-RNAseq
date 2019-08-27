#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=combine_GVCFs
#SBATCH --mem=40GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/scratch/Users/pein7187/Chickadee-RNAseq/new_pipeline/slurm_out_err/combine_GVCFs_all_%a.out
#SBATCH --error=/scratch/Users/pein7187/Chickadee-RNAseq/new_pipeline/slurm_out_err/combine_GVCFs_all_%a.err

module load gatk
module load samtools

#set -exuo pipefail

if [ $# -lt 1 ]
then
    echo "
    Combines and genotypes GVCFs, then filters the resulting VCF table, and removes indels.
    Input is all GVCF files generated from HaplotypeCaller.
    Output is final, filtered VCF table which will be used for downstream analyses.

    [-g] Path to directory containing GVCFs
    [-r] Path to ref genome fasta
    [-o] Path to output dir
    [-p] Output VCF prefix

    "

else
    while getopts g:r:o:p: option
    do
    case "${option}"
    in
    g) GVCFs_path=${OPTARG};;
    r) ref=${OPTARG};;
    o) outputdir=${OPTARG};;
    p) output_vcf_prefix=${OPTARG};;
    esac
    done

    for gvcf in $( ls $GVCFs_path*.gvcf ); do
    
        #if [ ! -f $vcf ];then
        #    echo "$vcf not found"
        #    exit
        #fi
    
        gvcf_args+="--variant $gvcf "
    
    done
    
    #combine gvcfs
    
    gatk CombineGVCFs \
    -R $ref \
    -O $outputdir"$output_vcf_prefix".joint.vcf \
    $gvcf_args
    
    
    #genotype gvcfs
    
    gatk GenotypeGVCFs \
    -R $ref \
    -O $outputdir"$output_vcf_prefix".vcf \
    -V $outputdir"$output_vcf_prefix".joint.vcf
    
    
   # #filter
   # 
   # gatk VariantFiltration \
   # -R $ref \
   # -V $outputdir"$output_vcf_prefix".vcf \
   # --window 35 \
   # --cluster 3 \
   # --filter-name "FS" \
   # --filter "FS > 30.0" \
   # --filter-name "QD" \
   # --filter "QD < 2.0" \
   # -O $outputdir"$output_vcf_prefix"_gatkfiltered.vcf
   # 

   # #remove indels

   # gatk SelectVariants \
   # -R $ref \
   # --variant $outputdir"$output_vcf_prefix"_gatkfiltered.vcf \
   # -O $outputdir"$output_vcf_prefix"_gatkfiltered_no_indels.vcf \
   # --select-type-to-exclude INDEL
   # 
   # #filter missing genotypes
   # #vcftools --max-missing 0.25 --gzvcf $OUT_PREFIX.vcf.gz --recode --stdout > $OUT_PREFIX\.vcf.gz

   # #remove intermediate files
   rm $outputdir"$output_vcf_prefix".joint.vcf
   # #rm $outputdir"$output_vcf_prefix".vcf
   # #rm $outputdir"$output_vcf_prefix"_gatkfiltered.vcf
fi
