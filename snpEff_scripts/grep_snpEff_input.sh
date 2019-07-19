#!/bin/bash

#SBATCH --job-name=grep
#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/scratch/Users/pein7187/slurm_out_err/moch_bcch_grep_snpEff_%a.out
#SBATCH --error=/scratch/Users/pein7187/slurm_out_err/moch_bcch_grep_snpEff_%a.err

#SNPs=/scratch/Users/pein7187/snpEff_analyses/BCCH/cach_bcch_allSNPs_no_nan.txt
SNPs2=/scratch/Users/pein7187/snpEff_analyses/BCCH/moch_bcch_allSNPs_no_nan.txt
vcf=/scratch/Users/pein7187/vcf/BCCH/joint_calling/all_chickadee_spp_BCCH_gatkfiltered.vcf.gz
outdir=/scratch/Users/pein7187/snpEff_analyses/BCCH/

#zgrep -w -f $SNPs $vcf > ${outdir}allSNPs_cach_bcch.snpEff_input.vcf
zgrep -w -f $SNPs2 $vcf > ${outdir}allSNPs_moch_bcch.snpEff_input.vcf
