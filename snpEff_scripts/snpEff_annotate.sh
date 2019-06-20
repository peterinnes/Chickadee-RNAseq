#!/bin/bash
path_to_snpEff=/home/peterinnes/snpEff/snpEff.jar
#outdir=/home/peterinnes/Chickadees/snpEff_analyses/BCCH/cach_vs_bcch/
#snpEff_db=taeGut_bcch_ann #database used for BCCH
snpEff_db=taeGut3.2.4.86 #database used for tG (zebrafinch)
#must be run in directory containing the snpEff input vcf files
for snpEff_input in *snpEff_input.vcf; do

prefix=$(echo $snpEff_input | sed -e 's/.snpEff_input.vcf//g')
java -jar $path_to_snpEff -v -c /home/peterinnes/snpEff/snpEff.config -stats $prefix $snpEff_db $snpEff_input > ${prefix}.snpEff_ann.vcf

mv $prefix ${prefix}_snpEff_summary.html
done
