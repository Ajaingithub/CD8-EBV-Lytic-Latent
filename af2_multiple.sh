#!/bin/bash
#SBATCH --job-name=test_multiple_job     # Job name
#SBATCH --gres=gpu:1                     # Request 1 GPU
#SBATCH --time=02:00:00                  # Max runtime: 2 hrs (adjust if needed)
#SBATCH --mem=32G                        # Memory request
#SBATCH --output=test_multiple_output_%A_%a.log  # Output log file
#SBATCH --error=gpu_step2_error_%A_%a.log  # Error log file

# Directories
DATA_DIR=/diazlab/data3/abhinav/resource/alphafold_database/
FASTA_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/fasta/fasta2/
OUTPUT_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/af2_output_structure/

# Activate conda environment
source ~/.bashrc
source /diazlab/data3/.abhinav/tools/miniconda3/etc/profile.d/conda.sh
conda activate alphafold2

# Get the list of FASTA files
FASTA_FILES=($(ls $FASTA_DIR))

# Determine how many FASTA files are there
NUM_FASTA_FILES=${#FASTA_FILES[@]}

# Ensure the SLURM array index is within bounds
if [ $SLURM_ARRAY_TASK_ID -ge $NUM_FASTA_FILES ]; then
    echo "SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) is out of bounds. Max index: $(($NUM_FASTA_FILES - 1))"
    exit 1
fi

# Select the FASTA file for this job
FASTA_FILE=${FASTA_FILES[$SLURM_ARRAY_TASK_ID]}

echo "Running for FASTA file: $FASTA_FILE"

# Run the AlphaFold job for the selected FASTA file
cd /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/alphafold-2.3.1/
bash run_alphafold.sh \
  -d $DATA_DIR \
  -o $OUTPUT_DIR \
  -f $FASTA_FILE \
  -t 2023-05-16

