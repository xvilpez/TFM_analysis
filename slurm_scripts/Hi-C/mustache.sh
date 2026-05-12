#!/bin/bash
#SBATCH --job-name=mustache
#SBATCH --mem=150gb
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=mustache_%A-%a.log

# ---------- Clean shell ----------
module purge
unset PYTHONPATH
unset LD_LIBRARY_PATH

# ---------- Activate Conda correctly ----------
module load Miniconda3/24.7.1-0
source /opt/ohpc/pub/libs/easybuild/4.5.0/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate mustache

# ---------- Sanity check ----------
which python
which mustache

# ---------- Variables ----------
describer=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
source ./config.sh

mkdir -p \
  "${restsite_folder}" \
  "${path_hicMatrix}" \
  "${path_coolMatrix}" \
  "${path_cooltools}" \
  "${path_loops}"

# ---------- Loop calling ----------
for number in 5 10 20
do
  echo "START mustache ${describer}"

  mustache \
    -f ${path_coolMatrix}/${describer}_${number}kb_KR.cool \
    -r ${number}kb -pt 0.1 \
    -o ${path_loops}/hic_${describer}_loops_01_${number}kb.tsv

  mustache \
    -f ${path_coolMatrix}/${describer}_${number}kb_KR.cool \
    -r ${number}kb -pt 0.05 \
    -o ${path_loops}/hic_${describer}_loops_05_${number}kb.tsv

  echo "END mustache ${describer}"
done
