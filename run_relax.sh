#!/bin/bash
#PBS -P va71
#PBS -q gpuvolta
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l mem=64GB
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/va71+scratch/va71
#PBS -l wd

# Load required modules
module load cuda/11.2.2

# Activate the AlphaFold environment
source /g/data/va71/Kerry/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/va71/Kerry/miniconda3/envs/alphafold

# Define paths
OUTPUT_DIR="/g/data/va71/Kerry/alphafold-2.3.1/output_test1"
PDB_PATH="$OUTPUT_DIR/unrelaxed_model_1_pred_0.pdb"
OUTPUT_PDB="$OUTPUT_DIR/relaxed_model_1.pdb"
LOG_FILE="$OUTPUT_DIR/relax.log"

# Create an empty log file
touch $LOG_FILE

# Run the relaxation step
python /g/data/va71/Kerry/alphafold-2.3.1/alphafold/relax/relax.py \
  --pdb_path=$PDB_PATH \
  --output_path=$OUTPUT_PDB \
  --max_iterations=0 \
  --use_gpu=True \
  &> $LOG_FILE

echo "Relaxation step completed. Check log at $LOG_FILE"