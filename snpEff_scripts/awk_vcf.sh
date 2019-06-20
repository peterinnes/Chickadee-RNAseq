#!/bin/bash
SNP_list=/home/peterinnes/Chickadees/snpEff_analyses/zebra_finch/cach_vs_bcch/cach_bcch_tG_fixdif.txt
vcf=~/Chickadees/vcf/all_chickadee_spp_tG_gatkfiltered.vcf.gz
prefix=fixdif_cach_bcch_tG
outdir=~/Chickadees/snpEff_analyses/zebra_finch/cach_vs_bcch/

#this script uses awk to match loci in one file (SNP_list) to a second file (vcf) in order to generate the requisite input file for snpEff annotations
awk 'BEGIN {OFS="\t"}NR==FNR{loci[$1$2];next}{if($1$2 in loci)print $1,$2,$3,$4,$5,$6,$7,$8}' $SNP_list \
<(gzip -dc $vcf) \
> ${outdir}${prefix}.snpEff_input.vcf
