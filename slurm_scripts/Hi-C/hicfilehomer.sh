#!/bin/sh

#SBATCH --job-name=homer
#SBATCH --mem=60gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=hicFile_%A-%a.txt
#SBATCH --array=1-6

# ========== VARIABLES ==========
names=( 697 HALO1 HiC-REH-EP1-1 HiC-REH-EP1-2 HiC-REH-EP1-Aro-1 RS411 )

describer=${names[$SLURM_ARRAY_TASK_ID-1]}
source ./config.sh

for dir in "${path_homer}" ; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done

# ========== MODULES ==========

module load homer/0.1
module load Java/17.0.2


echo "................................................................ START hic file ${describer} ................................................................"

tagDir2hicFile.pl ${path_homer}/${describer}_filtered  -juicer auto -genome ${genome} -juicerExe "java -jar /ijc/LABS/STIK/RAW/xavi/juicer_tools.1.9.9_jcuda.0.8.jar" -p 8

echo "................................................................ END hic file ${describer} ................................................................"

