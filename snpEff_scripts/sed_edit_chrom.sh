#!bin/bash/
input=/home/peterinnes/Chickadees/snpEff_analyses/BCCH/cach_vs_bcch/allSNPs_cach_bcch_BCCH.snpEff_input.vcf
outdir=/home/peterinnes/Chickadees/snpEff_analyses/BCCH/cach_vs_bcch/
prefix=allSNPs_cach_bcch_BCCH

cat $input \
| sed -e 's/^.*\(ref\)//g' -e 's/|//g' \
> ${outdir}${prefix}_chrom_edited.snpEff_input.vcf
