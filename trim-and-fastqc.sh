module load fastqc

samplesheet="/scratch/Users/pein7187/chickadee_fastq/chickadee_fastq_list_rows.txt"
fastq=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`

fastqc -t 16 $fastq --outdir /scratch/Users/pein7187/chickadee_fastq/fastqc

