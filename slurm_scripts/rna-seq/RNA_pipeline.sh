#!/bin/bash

#SBATCH --job-name=RNA-seq
#SBATCH --mem=50gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=RNA-seq_%A-%a.log
#SBATCH --array=1-12

# ========== VARIABLES ==========
names=(
  RNA-REH-EP1-Arg-1
RNA-REH-EP1-Arg-2
RNA-REH-EP1-Arg-3
RNA-REH-EP1-d189-1
RNA-REH-EP1-d189-2
RNA-REH-EP1-d189-3
RNA-REH-EP1-d89-1
RNA-REH-EP1-d89-2
RNA-REH-EP1-d89-3
RNA-REH-GFP-1
RNA-REH-GFP-2
RNA-REH-GFP-3
)

describer=${names[$SLURM_ARRAY_TASK_ID-1]}

path_fq='/ijc/LABS/STIK/RAW/rna_seq/rnaseq_marta_REHArg_mar26/fastq_files'
path_bam="/ijc/LABS/STIK/RAW/rna_seq/rnaseq_marta_REHArg_mar26/STAR"
path_bw="/ijc/LABS/STIK/RAW/rna_seq/rnaseq_marta_REHArg_mar26/bigwig"
genomeIndex='/mnt/beegfs/xvilchez/indexes/STAR'
GTFfile='/mnt/beegfs/public/references/annotation/human/gencode.v39.annotation.gtf'
effective_genome_size=2913022398

for dir in "${path_bam}" "${path_bw}" ; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done


# ========== MODULES ==========
module load fastqc-0.11.9-gcc-11.2.0-dd2vd2m
module load STAR/2.7.6a-GCC-11.2.0
module load samtools-1.12-gcc-11.2.0-n7fo7p2
module load deepTools/3.5.1-foss-2021b

# ========== STEP 1: FASTQC ==========

echo "................................................................ 1. START_FASTQC ${describer} ................................................................"

fastqc ${path_fq}/${describer}*.f*q.gz -o ${path_fq}

echo "................................................................ 1. END_FASTQC ${describer} ................................................................"

# ========== STEP 2: ALIGNMENT ==========

echo "................................................................ 3. START_ALIGNMENT ${describer} ................................................................"

STAR --genomeDir ${genomeIndex} \
    --sjdbGTFfile ${GTFfile} \
    --readFilesCommand zcat  \
    --readFilesIn ${path_fq}/${describer}_R1_*.fastq.gz ${path_fq}/${describer}_R2_*.fastq.gz \
    --runThreadN 4 \
    --outFileNamePrefix ${path_bam}/STAR_${describer} \
    --outSAMtype BAM SortedByCoordinate --outTmpDir tmp_${describer} \
    --outWigStrand Unstranded \
    --quantMode GeneCounts

samtools index ${path_bam}/STAR_${describer}Aligned.sortedByCoord.out.bam

echo "................................................................ 2. END_ALIGNMENT ${describer} ................................................................"


# ========== STEP 3: BAM TO BIGWIG =========

echo "................................................................ 9. START_BAMCOVERAGE ${describer} ................................................................"

bamCoverage --bam ${path_bam}/STAR_${describer}Aligned.sortedByCoord.out.bam \
        --outFileName ${path_bw}/${describer}.bw \
        --effectiveGenomeSize 2913022398 \
        --outFileFormat bigwig --binSize 1 \
        --normalizeUsing RPGC

echo "................................................................ 9. END_BAMCOVERAGE ${describer} ................................................................"
