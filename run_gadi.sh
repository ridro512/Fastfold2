#!/bin/bash
#PBS -P va71
#PBS -q gpuvolta
#PBS -l ngpus=1
#PBS -l ncpus=12
#PBS -l mem=128GB
#PBS -l walltime=48:00:00
#PBS -l storage=scratch/va71+gdata/va71
#PBS -l wd
#PBS -l jobfs=16GB

# Module loading
module load cuda/11.7.0

# -----------------------------------------------
# Setup Conda inside the PBS job
# -----------------------------------------------
source /g/data/va71/Kerry/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/va71/Kerry/miniconda3/envs/alphafold

# -----------------------------------------------
bash /g/data/va71/Kerry/alphafold-2.3.1/run_alphafold.sh \
    -d /g/data/va71/gaetan_data/alphafold_data/ \
    -o /g/data/va71/Kerry/alphafold-2.3.1/output_test1/ \
    -f /g/data/va71/Kerry/alphafold-2.3.1/input_test1/query.fasta \
    -t 2020-05-14 \
    -r false
