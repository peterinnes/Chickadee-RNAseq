#!/bin/bash

set pipefail -exuo

if [ $# -lt 1 ]
then
    echo "
    Post-VCF data processing and analysis. Takes VCF table through Fst calculations,
    SNP annotations w/ snpEff, and characterization of number of fixed differences in candidate pathways. All analysis is done for both pairwise species comparison: black-capped vs. mountain, black-capped vs. Carolina.
    
    [-v] Full path to gzipped VCF table containing all three chickadee species
    [-e] Full path to snpEff
    [-o] Full path to output directory
    [-s] Full path to samplesheet

    "
else
    while getopts v:e:o:s: option
    do
    case "${option}"
    in
    v) VCF=${OPTARG};;
    e) snpEff=${OPTARG};;
    o) outputdir={OPTARG};;
    s) samplesheet={OPTARG};;
    esac
    done

    if [ ! -d $outputdir ]; then
        mkdir $outputdir
        mkdir "$outputdir"moch_vs_bcch
        mkdir "$outputdir"cach_vs_bcch
    fi
        
    bcch="$outputdir"bcch_samples.txt
    moch="$outputdir"moch_samples.txt
    cach="$outputdir"cach_samples.txt
    
    #making population sample lists for vcftools Fst calculations b/w species pairs
    if [ ! -e $bcch ]; then
        cat $samplesheet | grep 'CU' > $bcch
    fi
    
    if [ ! -e $moch]; then
        cat $samplesheet | grep 'NK' > $moch
    fi
    
    if [ ! -e $cach]; then
        cat $samplesheet | grep 'B' > $cach
    fi
    
    
    #fst calculations
    vcftools --gzvcf $VCF --weir-fst-pop $moch --weir-fst-pop $bcch --stdout \
    > "$outputdir"moch_vs_bcch/moch_bcch.weir.fst
    vcftools --gzvcf $VCF --weir-fst-pop $cach --weir-fst-pop $bcch --stdout \
    > "$outputdir"cach_vs_bcch/cach_bcch.weir.fst
    
    #filter vcftools Fst output for fixed differences, then match these back to VCF to make snpEff_inputVCF
    #the sed command at the end is to make the chrom name of our vcf match the snpEff database chroms, 
    #since we are using annotations from an older zebra finch ref
    awk '$3==1' $outputdir/moch_vs_bcch/moch_bcch.weir.fst | awk 'BEGIN {OFS="\t"};{print$1,$2}' \
    | awk 'BEGIN {OFS="\t"}NR==FNR{loci[$1$2];next}{if($1$2 in loci)print $1,$2,$3,$4,$5,$6,$7,$8}' - <(gzip -dc $VCF) \
    | sed -e 's/^.*\(ref\)//g' -e 's/|//g' -e '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' \ 
    > "$outputdir"moch_vs_bcch/moch_bcch_fixdif.snpEff_input.vcf
    
    #same for cach_vs_bcch
    awk '$3==1' $outputdir/cach_vs_bcch/cach_bcch.weir.fst | awk 'BEGIN {OFS="\t"};{print$1,$2}' \
    | awk 'BEGIN {OFS="\t"}NR==FNR{loci[$1$2];next}{if($1$2 in loci)print $1,$2,$3,$4,$5,$6,$7,$8}' - <(gzip -dc $VCF) \
    | sed -e 's/^.*\(ref\)//g' -e 's/|//g' -e '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' \
    > "$outputdir"cach_vs_bcch/cach_bcch_fixdif.snpEff_input.vcf
    
   
    
    #snpEff
    java -jar "$snpEff"snpEff.jar -v -c "$snpEff"snpEff.config -stats moch_bcch taeGut_bcch_ann moch_bcch_fixdif_chrom_edited.snpEff_input.vcf \
    > moch_bcch_fixdif.ann.vcf
    
    java -jar "$snpEff"snpEff.jar -v -c "$snpEff"snpEff.config -stats cach_bcch taeGut_bcch_ann cach_bcch_fixdif_chrom_edited.snpEff_input.vcf \
    > cach_bcch_fixdif.snpEff_ann.vcf

