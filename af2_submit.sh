#!/bin/bash
#SBATCH --job-name=gpu_test_job    # Job name
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --time=02:00:00            # Max runtime: 2 hrs
#SBATCH --mem=16G                   # Memory request
#SBATCH --output=af2_CSASISGGDPYEQYF_output.log  # Output log file
#SBATCH --error=af2_CSASISGGDPYEQYF_error.log  # Error log file

# Activate your Conda environment
source ~/.bashrc
source /diazlab/data3/.abhinav/tools/miniconda3/etc/profile.d/conda.sh
source activate alphafold2
conda init
conda activate alphafold2    

cd /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/new_af/alphafold-2.3.1/
bash run_alphafold.sh \
  -d /diazlab/data3/abhinav/resource/alphafold_database/ \
    -o /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/af2_output_structure/ \
      -f  /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/fasta/CSASISGGDPYEQYF_1.fasta \
        -t 2020-05-16
