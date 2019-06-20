#Filter snpEff output for valid genes. First, filter out all loci out of chromosome range. Second, print $8, the annotation field of the vcf table. Third, find the gene name, which comes after "MODIFIER".
#!/bin/bash
#snpEff_output=/home/peterinnes/Chickadees/snpEff_analyses/BCCH/moch_vs_bcch/fixdif_moch_bcch.snpEff_ann.vcf
prefix=genes_fixdif_moch_bcch
#outdir=/home/peterinnes/Chickadees/snpEff_analyses/BCCH/moch_vs_bcch/

#must be run in the directory with snpEff output
for snpEff_output in *snpEff_ann*; do

prefix=$(echo $snpEff_output | sed -e 's/.snpEff_ann.vcf//' -e 's/^/genes_/')

cat $snpEff_output \
| grep -v '##\|OUT_OF_CHROMOSOME_RANGE' \
| cut -f 8 \
| sed -e 's/^.*\(MODIFIER\)//g' \
| awk 'BEGIN {FS="|"};{OFS="\t"};{print $2}' \
| sort -u \
| grep -v 'LOC\|id\|CHR_START' \
| tail -n +2 \
> ${prefix}.txt

done
