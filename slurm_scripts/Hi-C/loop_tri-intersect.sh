#!/bin/bash
#SBATCH --job-name=loop_tri-inter
#SBATCH --mem=15gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=loop_tri-intersect_%A-%a.txt
#SBATCH --array=1-1
module load bedtools2-2.30.0-gcc-11.2.0-q7z4zez

path="/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/loops"
mkdir -p ${path}/intersected_loops
cd ${path}

# # Prepare expanded bed files
# for describer in hic_HiC-REH-EP1 hic_HiC-REH-EP1-Aro hic_REH
# do
#     for number in 05_10kb
#     do
#         # expand loops and remove the first line
#         sed '1d' ${describer}_loops_${number}.tsv > ${describer}_loops_${number}.bed
#         awk '{
#             OFS="\t";
#             $2 = $2 - 10000;
#             $3 = $3 + 10000;
#             $5 = $5 - 10000;
#             $6 = $6 + 10000;
#             if ($2 < 0) $2 = 0;
#             if ($5 < 0) $5 = 0;
#             print $1, $2, $3, $4, $5, $6, $7, $8;
#         }' ${describer}_loops_${number}.bed > ${describer}_loops_${number}_exp30kb.bed
#     done
# done

for number in 05_10kb
do
    EP1="intersected_loops/high_confidence_EP1_loops_coords.bedpe"
    ARO="intersected_loops/high_confidence_ARO_loops_coords.bedpe"
    REH="hic_REH_loops_${number}_exp30kb.bed"

    OUT="${path}/intersected_loops"

    # Shared loops across all three conditions (EP1 ∩ ARO ∩ REH)
    pairToPair -a "$EP1" -b "$ARO" -type both | cut -f 1-6 | \
    pairToPair -a stdin -b "$REH" -type both \
    > "${OUT}/shared_loops_all3_${number}_exp30kb.bedpe"

    # EP1 only
    pairToPair -a "$EP1" -b "$ARO" -type notboth | cut -f 1-6 | \
    pairToPair -a stdin -b "$REH" -type notboth \
    > "${OUT}/EP1_specific_only_${number}_exp30kb.bedpe"

    # ARO only
    pairToPair -a "$ARO" -b "$EP1" -type notboth | cut -f 1-6 | \
    pairToPair -a stdin -b "$REH" -type notboth \
    > "${OUT}/ARO_specific_only_${number}_exp30kb.bedpe"

    # REH only
    pairToPair -a "$REH" -b "$EP1" -type notboth | cut -f 1-6 | \
    pairToPair -a stdin -b "$ARO" -type notboth \
    > "${OUT}/REH_specific_only_${number}_exp30kb.bedpe"

    # Shared EP1 + ARO only
    pairToPair -a "$EP1" -b "$ARO" -type both | cut -f 1-6 | \
    pairToPair -a stdin -b "$REH" -type notboth \
    > "${OUT}/shared_EP1_ARO_only_${number}_exp30kb.bedpe"

    # Shared EP1 + REH only
    pairToPair -a "$EP1" -b "$REH" -type both | cut -f 1-6 | \
    pairToPair -a stdin -b "$ARO" -type notboth \
    > "${OUT}/shared_EP1_REH_only_${number}_exp30kb.bedpe"

    # Shared ARO + REH only
    pairToPair -a "$ARO" -b "$REH" -type both | cut -f 1-6 | \
    pairToPair -a stdin -b "$EP1" -type notboth \
    > "${OUT}/shared_ARO_REH_only_${number}_exp30kb.bedpe"
done
