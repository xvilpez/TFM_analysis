#!/bin/bash
#SBATCH --job-name=STAR_index_GRCh38
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00
#SBATCH --output=STAR_index_GRCh38.log

module load STAR/2.7.6a-GCC-11.2.0

STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir /mnt/beegfs/xvilchez/indexes/STAR \
  --genomeFastaFiles /mnt/beegfs/public/references/genome/human/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile /mnt/beegfs/public/references/annotation/human/gencode.v39.annotation.gtf \
  --sjdbOverhang 149
