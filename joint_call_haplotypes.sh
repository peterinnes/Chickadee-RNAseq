module load picard
module load gatk
module load samtools

samplesheet=/scratch/Users/pein7187/chickadee_samplesheet.txt
sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
bam_path=/scratch/Users/pein7187/sorted_RGadded_dupmarked_split_bams/zebra_finch/
bam_suffix=_trimmed_tG.sorted_RGadded_dupmarked_split.bam
ref=/scratch/Users/pein7187/taeGut/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa
outdir=/scratch/Users/pein7187/vcf/zebra_finch/joint_calling/

gatk HaplotypeCaller \
-R $ref \
-I $bam_path"$sample"$bam_suffix \
-O $outdir"$sample".gvcf \
-ERC GVCF \
--dont-use-soft-clipped-bases \
--standard-min-confidence-threshold-for-calling 20
