#  Folder paths
path_fq='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/fastq_files'
path_bam='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/bowtie2_results'
restsite_folder='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/rest_sites'
path_hicMatrix='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/hicMatrix'
path_homer='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/homer'
path_coolMatrix='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/coolMatrix'
path_loops='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/loops'
path_cooltools='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/cooltools'
path_cscore='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/cscore'

# Files and programs
HICUP_trunc='../bin/HiCUP-0.9.2/hicup_truncater'   # requirements
indexgenome='/mnt/beegfs/public/references/index/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as'   # requirements
refgenome='/mnt/beegfs/public/references/genome/human/GRCh38.primary_assembly.genome.fa'
ref_compartments='ref/refcool_MRC5_ATAC_100kb.bed'  # requirements
genome100kb='../files/hg38_100kb.bed'
genome_RabI='../files/hg38_100kb_RabI_clean.bed'
blacklist='../files/black_list'

# Parameters
N=$(wc -l < samples.txt)
restriction_enzyme='^GATC,MboI'
restrictionSequence='GATC'
genome='hg38'
