module load samtools
module load picard
module load gatk

samplesheet="/scratch/Users/pein7187/chickadee_samplesheet.txt"
sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
path_to_sams=/scratch/Users/pein7187/STAR_mapping_2passMode/trimmed/zebra_finch/sams/
bamoutdir=/scratch/Users/pein7187/sorted_RGadded_dupmarked_split_bams/zebra_finch/
ref=/scratch/Users/pein7187/taeGut/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa
ref_abrv=tG

samtools view -@ 12 -S -b $path_to_sams"$sample"_trimmedAligned.out.sam > $bamoutdir"$sample"_trimmed_$ref_abrv\.bam
samtools sort -@ 12 -T temp -o $bamoutdir"$sample"_trimmed_$ref_abrv\.sorted.bam $bamoutdir"$sample"_trimmed_$ref_abrv\.bam
rm $bamoutdir"$sample"_trimmed_$ref_abrv\.bam

java -jar /opt/picard/2.6.0/picard-2.6.0.jar AddOrReplaceReadGroups \
I=$bamoutdir"$sample"_trimmed_$ref_abrv\.sorted.bam \
O=$bamoutdir"$sample"_trimmed_$ref_abrv\.sorted_RGadded.bam \
RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample
rm $bamoutdir"$sample"_trimmed_$ref_abrv\.sorted.bam

java -jar /opt/picard/2.6.0/picard-2.6.0.jar MarkDuplicates \
I=$bamoutdir"$sample"_trimmed_$ref_abrv\.sorted_RGadded.bam \
O=$bamoutdir"$sample"_trimmed_$ref_abrv\.sorted_RGadded_dupmarked.bam \
M=$bamoutdir"$sample"_trimmed_$ref_abrv\.dupmarked_metrics.txt
rm $bamoutdir"$sample"_trimmed_$ref_abrv\.sorted_RGadded.bam

gatk SplitNCigarReads -R $ref -I $bamoutdir"$sample"_trimmed_$ref_abrv\.sorted_RGadded_dupmarked.bam -O $bamoutdir"$sample"_trimmed_$ref_abrv\.sorted_RGadded_dupmarked_split.bam
rm $bamoutdir"$sample"_trimmed_$ref_abrv\.sorted_RGadded_dupmarked.bam

