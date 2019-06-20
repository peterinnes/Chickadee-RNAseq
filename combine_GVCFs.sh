#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=combine_GVCFs
#SBATCH --mem=40GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/scratch/Users/pein7187/slurm_out_err/tG_combine_GVCFs_all_%a.out
#SBATCH --error=/scratch/Users/pein7187/slurm_out_err/tG_combine_GVCFs_all_%a.err

ref=/scratch/Users/pein7187/taeGut/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa
vcf_csv=/scratch/Users/pein7187/vcf/zebra_finch/joint_calling/all_spp_gvcf_list_tG.csv
outputdir=/scratch/Users/pein7187/vcf/zebra_finch/joint_calling/
output_vcf_prefix=all_chickadee_spp_tG

module load gatk
module load samtools

set -exuo pipefail

#echo "STARTING reference=$ref vcf_csv=$vcf_csv output_vcf_prefix=$output_vcf_prefix $0"

vcf_args=""
for vcf in $( cat $vcf_csv | tr ',' '\n' ); do

    if [ ! -f $vcf ];then
        echo "$vcf not found"
        exit
    fi

    vcf_args+="--variant $vcf "

done

gatk CombineGVCFs -R $ref -O $outputdir"$output_vcf_prefix".joint.vcf $vcf_args

gatk GenotypeGVCFs -R $ref -O $outputdir"$output_vcf_prefix".vcf -V $outputdir"$output_vcf_prefix".joint.vcf

gatk VariantFiltration \
-R $ref \
-V $outputdir"$output_vcf_prefix".vcf \
--window 35 \
--cluster 3 \
--filter-name "FS" \
--filter "FS > 30.0" \
--filter-name "QD" \
--filter "QD < 2.0" \
-O $outputdir"$output_vcf_prefix"_gatkfiltered.vcf

rm $outputdir"$output_vcf_prefix".joint.vcf
rm $outputdir"$output_vcf_prefix".vcf