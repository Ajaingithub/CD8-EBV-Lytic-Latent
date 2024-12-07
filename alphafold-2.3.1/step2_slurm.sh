#!/bin/bash
#SBATCH --job-name=gpu_test_job    # Job name
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --time=02:00:00            # Max runtime: 2 hrs (adjust if needed)
#SBATCH --mem=50G                   # Memory request
#SBATCH --output=gpu_step2_output.log  # Output log file
#SBATCH --error=gpu_step2_error.log  # Error log file

bash step2_2.sh
