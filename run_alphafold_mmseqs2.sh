#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l jobfs=400GB
#PBS -q normal
#PBS -P va71
#PBS -l storage=gdata/va71

# Load required modules
module load python3/3.9.2
module load cuda/11.2.2

# Set up environment
export CUDA_VISIBLE_DEVICES=0
export ALPHAFOLD_DATA_DIR=/g/data/va71/gaetan_data/alphafold_data

# Change to AlphaFold directory
cd /g/data/va71/Kerry/alphafold-2.3.1

# Activate virtual environment if you have one
# source venv/bin/activate

# Run AlphaFold with MMseqs2
python run_alphafold.py \
    --fasta_paths=$1 \
    --max_template_date=2023-05-14 \
    --model_preset=monomer \
    --db_preset=full_dbs \
    --data_dir=/g/data/va71/gaetan_data/alphafold_data \
    --output_dir=./output \
    --num_multimer_predictions_per_model=1
