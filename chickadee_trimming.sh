module load java
module load trimmomatic

samplesheet="/scratch/Users/pein7187/chickadee_fastq/chickadee_fastq_list_rows.txt"
fastq=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}' | sed -e 's/.fastq.gz//g' -e 's/^.\{45\}//'`

java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar SE -threads 16 -phred33 $fastq /scratch/Users/pein7187/chickadee_fastq/trimmed_reads/${sample}.trimmed.fastq.gz ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/TruSeq2-SE.fa:2:30:10 HEADCROP:15 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36
