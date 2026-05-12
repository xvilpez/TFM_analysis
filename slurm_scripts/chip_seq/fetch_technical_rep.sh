#!/bin/bash

#SBATCH --job-name=fq_dump
#SBATCH --mem=80gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=pd_%A-%a.txt

#MODULS NECESSARIS

module load SRA-Toolkit/2.10.9-gompi-2021b

mkdir -p ../data/raw

cd ../data/raw

for describer in SRR60963{0..1}; do

	prefetch --max-size 80GB -v ${describer}

	fastq-dump ${describer}/*.sra --gzip --outdir .

done

cat SRR609630.fastq.gz SRR609631.fastq.gz > CD19.RO_01736.fastq.gz

