#!/bin/bash

#SBATCH --job-name=quality_control
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH --output=QC_pretrimm%j.log

#required modules

module load fastqc-0.11.9-gcc-11.2.0-dd2vd2m

WD=/ijc/LABS/STIK/RAW/chip_seq/Xavi_TFM/REH_public_marta

cd $WD/fastq_files

for fastq in *.fastq.gz
do

echo "############################################################ START QUALITY CONTROL ${fastq} ###########################################################################"

fastqc -o . $fastq

echo "############################################################# END QUALITY CONTROL ${fastq} ###########################################################################"

done
