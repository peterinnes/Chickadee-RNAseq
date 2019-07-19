#used this script in the per-sample variant calling pipeline. produces vcfs much smaller in size than the joint-calling pipeline.

#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=CombineVariants
#SBATCH --mem=40GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/scratch/Users/pein7187/slurm_out_err/combine_variants_%a.out
#SBATCH --error=/scratch/Users/pein7187/slurm_out_err/combine_variants_%a.err

module load gatk/3.7.0 #using GATK 3.7 because CombineVariants is not in GATK 4, and MergeVCFs can only be used to combine VCFs from the same sample.
#however, GATK discourages merging single sample VCFs (per Sheila on GATK forum)--furthermore the final VCF resulting from this method resulted in way more missing genotype data than "joint-calling," so I did that instead; see combine_GVCFs script.

#moch_bcch_vcf_list=/scratch/Users/pein7187/vcf/BCCH/moch_bcch_vcf.list
#cach_bcch_vcf_list=/scratch/Users/pein7187/vcf/BCCH/cach_bcch_vcf.list
outdir=/scratch/Users/pein7187/vcf/BCCH/
ref=/scratch/Users/pein7187/BCCH_RefGenome/BCCH_pseudochr_rename.fasta

#cach_bcch
java -jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $ref \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_B_77292_GTCCGC_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_B_77293_GTGAAA_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_B_77294_CCGTCC_L004_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_B_77295_GTCCGC_L004_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_B_77296_GTGAAA_L004_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_CU54780_CGATGT_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_CU54781_TGACCA_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_CU54782_ACAGTG_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_CU54794_GCCAAT_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_CU54795_CAGATC_L003_R1_001_filtered.vcf \
-o ${outdir}cach_bcch_combined.vcf \
-genotypeMergeOptions UNIQUIFY

#moch_bcch
java -jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $ref \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_CU54780_CGATGT_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_CU54781_TGACCA_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_CU54782_ACAGTG_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_CU54794_GCCAAT_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_CU54795_CAGATC_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_NK218349_CTTGTA_L004_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_NK218350_AGTCAA_L004_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_NK218351_AGTTCC_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_NK218352_ATGTCA_L003_R1_001_filtered.vcf \
-V /scratch/Users/pein7187/vcf/BCCH/Chickadee_NK218353_CCGTCC_L003_R1_001_filtered.vcf \
-o ${outdir}moch_bcch_combined.vcf \
-genotypeMergeOptions UNIQUIFY
