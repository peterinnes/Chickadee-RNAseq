#!/bin/bash

tG_vcf=/scratch/Users/pein7187/Chickadee-RNAseq/vcf/zebra_finch/all_chickadee_spp_tG_gatkfiltered.vcf.gz
BCCH_vcf=/scratch/Users/pein7187/Chickadee-RNAseq/vcf/BCCH/joint_calling/all_chickadee_spp_BCCH_gatkfiltered.vcf.gz
bcch=/scratch/Users/pein7187/Chickadee-RNAseq/bcch_samples.txt
moch=/scratch/Users/pein7187/Chickadee-RNAseq/moch_samples.txt
cach=/scratch/Users/pein7187/Chickadee-RNAseq/cach_samples.txt
outdir=/scratch/Users/pein7187/Chickadee-RNAseq/vcf/BCCH/joint_calling/vcftools_analysis/

module load vcftools

vcftools --gzvcf $BCCH_vcf --weir-fst-pop $moch --weir-fst-pop $bcch --fst-window-size 100000 --fst-window-step 100000 --stdout > $outdir\sliding_window100000_moch_bcch_BCCH.weir.fst 
vcftools --gzvcf $BCCH_vcf --weir-fst-pop $cach --weir-fst-pop $bcch --fst-window-size 100000 --fst-window-step 100000 --stdout > $outdir\sliding_window100000_cach_bcch_BCCH.weir.fst 
