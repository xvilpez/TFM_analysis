#!/bin/bash
 
#SBATCH --job-name=fq_dump
#SBATCH --mem=80gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=pd_%A-%a.txt
#SBATCH --array=1-1

#MODULS NECESSARIS

module load SRA-Toolkit/2.10.9-gompi-2021b

#mkdir data
#mkdir fastq_files

names=( SRR8550460 )

describer=${names[$SLURM_ARRAY_TASK_ID-1]}

cd ../data/raw 

prefetch --max-size 80GB -v ${describer} 

fastq-dump ${describer}/*.sra --gzip --outdir .
