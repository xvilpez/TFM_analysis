#!/bin/bash

#SBATCH --job-name=heatmap
#SBATCH --mem=5gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=6
#SBATCH --output=heatmap_%A-%a.log

module load deepTools/3.5.1-foss-2021b

path=/ijc/LABS/STIK/RAW/chip_seq/Xavi_TFM

cd $path

mkdir -p heatmaps

# computeMatrix reference-point \
#     -S bigwig/E-MTAB-8177:WT-RUNX1.bw bigwig/ChIP-seq-NALM6-a-RUNX1-1.bw bigwig/ChIP-seq-NALM6-a-RUNX1-2.bw bigwig/Nalm6_Runx1_rep1.bw bigwig/Nalm6_Runx1_rep2.bw bigwig/REH_RUNX1_rep1.bw bigwig/REH_RUNX1_rep2.bw \
#     -R 697HF_E2APBX1-HF_peaks.bed \
#     -a 5000 -b 5000 \
#     -p 6 \
#     --missingDataAsZero \
#     --referencePoint center \
#     --skipZeros \
#     -out heatmaps/RUNX1.mat.gz

# plotHeatmap -m heatmaps/RUNX1.mat.gz \
#     --colorList 'white,blue' 'white,red' 'white,red' 'white,red' 'white,red' 'white,green' 'white,green' \
#     -x "Peak center distance (bp)" \
#     --regionsLabel "TCF3::PBX1 peak regions" \
#     --samplesLabel "697" "Nalm6-1" "Nalm6-2" "Nalm6_rep1" "Nalm6_rep2" "REH_rep1" "REH_rep2" \
#     --plotTitle "RUNX1 ChIP-seq Heatmap" \
#     -out heatmaps/RUNX1.pdf

# computeMatrix reference-point \
#     -S bigwig/H3K27Ac_ChIP_sample_9.bw bigwig/H3K27Ac_ChIP_sample_8.bw bigwig/H3K27Ac_ChIP_sample_7.bw bigwig/H3K27Ac_ChIP_sample_1.bw bigwig/H3K27Ac_ChIP_sample_2.bw bigwig/H3K27Ac_ChIP_sample_3.bw bigwig/H3K27Ac_ChIP_sample_4.bw bigwig/H3K27Ac_ChIP_sample_5.bw bigwig/H3K27Ac_ChIP_sample_6.bw  bigwig/CD19.RO_01736.bw \
#     -R 697HF_E2APBX1-HF_peaks.bed \
#     -a 5000 -b 5000 \
#     -p 6 \
#     --missingDataAsZero \
#     --referencePoint center \
#     --skipZeros \
#     -o heatmaps/H3K27Ac_samples.mat.gz

plotHeatmap -m heatmaps/H3K27Ac_samples.mat.gz \
    --colorList 'white,blue' 'white, red' 'white, green' 'white, orange' 'white, orange' 'white, purple' 'white, purple' 'white, purple' 'white, purple' 'white, lightskyblue'\
    -x "Peak center distance (bp)" \
    --regionsLabel "TCF3::PBX1 peak regions" \
    --samplesLabel "TCF3::PBX1" "Ph_like" "PAX5" "BCR::ABL1" "BCR::ABL1" "KMT2A" "KMT2A" "KMT2A" "KMT2A" "CD19+ Control" \
    --plotTitle "H3K27Ac Patients ChIP-seq Heatmap" \
    -out heatmaps/H3K27Ac_samples.pdf \

# computeMatrix reference-point \
#     -S bigwig/697_H3K27ac_ChIP-seq.bw bigwig/Nalm6_H3K27ac.bw bigwig/ChIP-seq-NALM6-a-H3K27ac-1.bw bigwig/ChIP-seq-NALM6-a-H3K27ac-2.bw bigwig/REH_H3K27ac.bw  bigwig/REH_H3K27ac_ChIP-seq.bw bigwig/H000H3K27acX*.bw bigwig/H168H3K27acX*.bw bigwig/CD19.RO_01736.bw \
#     -R 697HF_E2APBX1-HF_peaks.bed \
#     -a 5000 -b 5000 \
#     -p 6 \
#     --missingDataAsZero \
#     --referencePoint center \
#     --skipZeros \
#     -o heatmaps/H3K27Ac_cells.mat.gz

# plotHeatmap -m heatmaps/H3K27Ac_cells.mat.gz \
#     --colorList 'white,blue' 'white, red' 'white, red' 'white, red' 'white, green' 'white, green' 'white, purple' 'white, purple' 'white, purple' 'white, purple' 'white, lightskyblue'\
#     -x "Peak center distance (bp)" \
#     --regionsLabel "TCF3::PBX1 peak regions" \
#     --samplesLabel "697" "Nalm6" "NALM6-1" "NALM6-2" "REH" "REH_ChIP-seq" "RCH-1" "RCH-2" "iMac-1" "iMac-2" "CD19+_Bcell" \
#     --plotTitle "H3K27Ac Cells ChIP-seq Heatmap" \
#     -out heatmaps/H3K27Ac_cells.pdf \
