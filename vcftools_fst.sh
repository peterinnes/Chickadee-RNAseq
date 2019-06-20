set pipefail -exuo

tG_vcf=/home/peterinnes/Chickadees/vcf/all_chickadee_spp_tG_gatkfiltered.vcf.gz
BCCH_vcf=/home/peterinnes/Chickadees/vcf/all_chickadee_spp_BCCH_gatkfiltered.vcf.gz
bcch=/home/peterinnes/Chickadees/bcch_samples.txt
moch=/home/peterinnes/Chickadees/moch_samples.txt
cach=/home/peterinnes/Chickadees/cach_samples.txt
outdir=/home/peterinnes/Chickadees

#making population sample lists
if [ ! -e $bcch ]; then
    cat /scratch/Users/pein7187/chickadee_samplesheet.txt | grep 'CU' > $bcch
fi

if [ ! -e $moch]; then
    cat /scratch/Users/pein7187/chickadee_samplesheet.txt | grep 'NK' > $moch
fi

if [ ! -e $cach]; then
    cat /scratch/Users/pein7187/chickadee_samplesheet.txt | grep 'B' > $cach
fi

if [ ! -d $outdir ]; then
    mkdir $outdir
fi


#fst calculations for data aligned to zebra finch genome
vcftools --gzvcf $tG_vcf --weir-fst-pop $moch --weir-fst-pop $bcch --stdout > $outdir/moch_bcch_tG.weir.fst
vcftools --gzvcf $tG_vcf --weir-fst-pop $cach --weir-fst-pop $bcch --stdout > $outdir/cach_bcch_tG.weir.fst

#fst calculations for data aligned to black-capped chickadee genome
vcftools --gzvcf $BCCH_vcf --weir-fst-pop $moch --weir-fst-pop $bcch --stdout > $outdir/moch_bcch_BCCH.weir.fst
vcftools --gzvcf $BCCH_vcf --weir-fst-pop $cach --weir-fst-pop $bcch --stdout > $outdir/cach_bcch_BCCH.weir.fst
