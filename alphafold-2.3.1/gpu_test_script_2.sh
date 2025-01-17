#!/bin/bash
#SBATCH --job-name=gpu_test_job    # Job name
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --time=00:10:00            # Max runtime: 10 minutes (adjust if needed)
#SBATCH --mem=8G                   # Memory request
#SBATCH --output=gpu_test_output_2.log  # Output log file
#SBATCH --error=gpu_test_error_2.log  # Error log file

# Load necessary modules (make sure the paths are correct for your HPC)
# module load cuda/12.4              # Load CUDA module (adjust if needed)
# module load python/3.8             # Load Python module (adjust if needed)

# Activate your Conda environment
# Use your specific Conda environment
source ~/.bashrc
source /diazlab/data3/.abhinav/tools/miniconda3/etc/profile.d/conda.sh
source activate alphafold
conda init
conda activate alphafold

# Python script to test GPU functionality
python - <<EOF
import jax
import jax.numpy as jnp

# Check available devices
devices = jax.devices()
print("Available devices:", devices)

# Create a matrix of ones and perform a simple operation
x = jnp.ones((1000, 1000))  # Create a 1000x1000 matrix of ones
y = jnp.dot(x, x)           # Matrix multiplication (dot product)

# Print a small part of the result to confirm
print("Result of matrix multiplication (partial):", y[:5, :5]) # Print top-left 5x5 block
EOF

bash run_alphafold.sh -d /diazlab/data3/abhinav/resource/alphafold_database/ -o ./dummy_test_slurm_gpu/ -f /diazlab/data3/abhinav/resource/alphafold_database/SSYTGSRTLV_8844.fasta -t 2024-12-01

# End of script

# Check if GPU is available
if tf.test.is_gpu_available():
    print("GPU is available")
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        for gpu in gpus:
            print(f"GPU Name: {tf.config.experimental.get_device_details(gpu)['device_name']}")
            print(f"Total Memory: {tf.config.experimental.get_memory_info(gpu)['current'] / 1024**2} MB")
else:
    print("GPU is not available.")

import tensorflow as tf

# Check if GPU is available
if tf.test.is_gpu_available():
    print("GPU is available")
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        # Set memory growth to prevent TensorFlow from allocating all GPU memory
        tf.config.experimental.set_memory_growth(gpus[0], True)
        print(f"GPU Name: {tf.config.experimental.get_device_details(gpus[0])['device_name']}")
        print(f"Total Memory: {tf.config.experimental.get_memory_info(gpus[0])['current'] / 1024**2} MB")
else:
    print("GPU is not available.")

import tensorflow as tf

conda create --name alphafold_gpu python=3.8
conda activate alphafold_gpu
