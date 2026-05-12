#!/bin/bash
 
#SBATCH --job-name=fq_dump
#SBATCH --mem=80gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=pd_%A-%a.txt
#SBATCH --array=1-1

#MODULS NECESSARIS

module load SRA-Toolkit/2.10.9-gompi-2021b

WD=/ijc/LABS/STIK/RAW/chip_seq/Xavi_TFM/REH_public_marta

mkdir -p $WD/fastq_files
cd $WD/fastq_files

names=(
#SRR26065561
#SRR26065560
SRR26065556 
)

describer=${names[$SLURM_ARRAY_TASK_ID-1]} 

prefetch --max-size 80GB -v ${describer} 

fastq-dump --split-files ${describer}/*.sra --gzip --outdir .
