#per-sample
module load picard
module load gatk
module load samtools

samplesheet=/scratch/Users/pein7187/chickadee_samplesheet.txt
sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
bam_path=/scratch/Users/pein7187/sorted_RGadded_dupmarked_split_bams/BCCH/
bam_suffix=_trimmed_bcchAligned.sorted_RGadded_dupmarked_split.bam
ref=/scratch/Users/pein7187/BCCH_RefGenome/BCCH_pseudochr_rename.fasta
outdir=/scratch/Users/pein7187/vcf/BCCH/

gatk HaplotypeCaller \
-R $ref \
-I $bam_path"$sample"$bam_suffix \
-O $outdir"$sample".vcf \
-dont-use-soft-clipped-bases \
--standard-min-confidence-threshold-for-calling 20
#gzip $outdir"$sample".vcf
#tabix $outdir"$sample".vcf.gz

gatk VariantFiltration \
--R $ref \
--V $outdir"$sample".vcf \
--window 35 \
--cluster 3 \
--filter-name "FS" \
--filter "FS > 30.0" \
--filter-name "QD" \
--filter "QD < 2.0" \
-O $outdir"$sample"_filtered.vcf
#tabix $outdir"$sample"_filtered.vcf.gz


