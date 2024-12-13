# sinfo -o "%P %G %D %N %C" ## To check for the Cuda partition
# srun --partition=common --gres=gpu:1 --time=2:00:00 --pty bash
# srun --partition=common --gres=gpu:1 --time=2:00:00 --pty bash
# nvidia-smi
## installing https://github.com/TencentAILabHealthcare/DeepAIR/tree/main
# pip install tensorflow
# module load cuda/11.7

# (DeepAIR) [ajain@c4-n39:job=244467 spatial]$ nvcc --version
# nvcc: NVIDIA (R) Cuda compiler driver
# Copyright (c) 2005-2022 NVIDIA Corporation
# Built on Wed_Jun__8_16:49:14_PDT_2022
# Cuda compilation tools, release 11.7, V11.7.99
# Build cuda_11.7.r11.7/compiler.31442593_0

### Install pytorch
# pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu117

import torch
print(torch.__version__)  # Prints the installed PyTorch version
print(torch.cuda.is_available())  # Should return True if CUDA is accessible
print(torch.cuda.get_device_name(0))  # Prints the name of the GPU

import os
from transformers import AutoModel, AutoTokenizer

os.chdir("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent")
model_name = "Rostlab/prot_bert_bfd"
output_dir = "./ProtTrans/prot_bert_bfd"

# Download and save model and tokenizer
model = AutoModel.from_pretrained(model_name)
tokenizer = AutoTokenizer.from_pretrained(model_name)

# Save files to the specified directory
model.save_pretrained(output_dir)
tokenizer.save_pretrained(output_dir)

print(f"Model and tokenizer saved in: {output_dir}")

#### Running the DeepAIR
### Preprocessing the files
## Step 1
python step1.py \
    --AIR_file_path ./sampledata/CoV-AbDab_example.csv \ 
    --output_table ./sampledata/CoV-AbDab_example_CDR3Region.csv \ 
    --output_fasta_folder ./sampledata/fasta \ 

# step2.py
# export DBS=/diazlab/data3/abhinav/resource/alphafold_database
# singularity exec --writable-tmpfs \
# --bind /diazlab,/francislab,/scratch \
# ${DBS}/AlphaFold.sif \
# /app/run_alphafold.sh \
# --use_gpu_relax=False \
# --bfd_database_path=${DBS}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
# --uniref30_database_path=${DBS}/uniref30/UniRef30_2021_03 \
# --pdb70_database_path=${DBS}/pdb70/pdb70 \
# --uniref90_database_path=${DBS}/uniref90/uniref90.fasta \
# --mgnify_database_path=${DBS}/mgnify/mgy_clusters_2022_05.fa \
# --template_mmcif_dir=${DBS}/pdb_mmcif/mmcif_files/ \
# --obsolete_pdbs_path=${DBS}/pdb_mmcif/obsolete.dat \
# --data_dir=${DBS}/ \
# --max_template_date=3000-01-01 \
# --model_preset=monomer \
# --output_dir=${PWD}/out/ \
# --fasta_paths=SSYTGSRTLV_8844.fasta

## with using GPUs
export DBS=/diazlab/data3/abhinav/resource/alphafold_database
singularity exec --nv --writable-tmpfs \
--bind /diazlab,/francislab,/scratch \
${DBS}/AlphaFold.sif \
/app/run_alphafold.sh \
--use_gpu_relax=False \
--bfd_database_path=${DBS}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
--uniref30_database_path=${DBS}/uniref30/UniRef30_2021_03 \
--pdb70_database_path=${DBS}/pdb70/pdb70 \
--uniref90_database_path=${DBS}/uniref90/uniref90.fasta \
--mgnify_database_path=${DBS}/mgnify/mgy_clusters_2022_05.fa \
--template_mmcif_dir=${DBS}/pdb_mmcif/mmcif_files/ \
--obsolete_pdbs_path=${DBS}/pdb_mmcif/obsolete.dat \
--data_dir=${DBS}/ \
--max_template_date=3000-01-01 \
--model_preset=monomer \
--output_dir=${PWD}/gpu_out_2/ \
--fasta_paths=SSYTGSRTLV_8844.fasta

### Since we have to make some changes in the model.py. We cannot use the singularity image of the Jake. 
# I have used Kaliana lab to perform local alpha fold installation
# https://github.com/kalininalab/alphafold_non_docker
# cd /diazlab/data3/.abhinav/tools/miniconda3/envs/deepair/lib/python3.8/site-packages/ && patch -p0 < $alphafold_path/docker/openmm.patch
# # $alphafold_path variable is set to the alphafold git repo directory (absolute path)
# conda install -y -c conda-forge openmm==7.5.1 pdbfixer
# conda install nvidia/label/cuda-12.4.1::cuda-toolkit --no-channel-priority
# conda install -y -c bioconda hmmer hhsuite==3.3.0 kalign2
# pip install nvidia-cudnn-cu12==8.9.4.25
# pip install nvidia-cublas-cu12==12.4.5.8
# pip install --upgrade --no-cache-dir jax==0.4.11 jaxlib==0.4.11+cuda12.cudnn88 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
# pip install tensorflow[and-cuda]


### Test if GPU is being able to called by JAX and tensorflow
import jax
import jax.numpy as jnp

# Check available devices
devices = jax.devices()
print("Available devices:", devices)

# Create a matrix of ones and perform a simple operation
x = jnp.ones((1000, 1000))  # Create a 1000x1000 matrix of ones
y = jnp.dot(x, x)           # Matrix multiplication (dot product)

# Print a small part of the result to confirm
print("Result of matrix multiplication (partial):", y[:5, :5])  # Print top-left 5x5 block

# cd /diazlab/data3/.abhinav/tools/miniconda3/envs/alphafold2/lib/python3.8/site-packages/ && patch -p0 < $alphafold_path/docker/openmm.patch

### Checking if the GPU is working
#!/bin/bash
#SBATCH --job-name=gpu_test_job    # Job name
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --time=02:00:00            # Max runtime: 10 minutes (adjust if needed)
#SBATCH --mem=8G                   # Memory request
#SBATCH --output=gpu_test_output.log  # Output log file
#SBATCH --error=gpu_test_error.log  # Error log file

# Load necessary modules (make sure the paths are correct for your HPC)
# module load cuda/12.4              # Load CUDA module (adjust if needed)
# module load python/3.8             # Load Python module (adjust if needed)

# Activate your Conda environment
source activate alphafold           # Use your specific Conda environment

# Python script to test GPU functionality
python - << EOF
import jax
import jax.numpy as jnp

# Check available devices
devices = jax.devices()
print("Available devices:", devices)

# Create a matrix of ones and perform a simple operation
x = jnp.ones((1000, 1000))  # Create a 1000x1000 matrix of ones
y = jnp.dot(x, x)           # Matrix multiplication (dot product)

# Print a small part of the result to confirm
print("Result of matrix multiplication (partial):", y[:5, :5])  # Print top-left 5x5 block
EOF

# End of script
### Running Step2
DATA_DIR=/diazlab/data3/abhinav/resource/alphafold_database/   # Path to the directory containing the AlphaFold 2 downloaded data.
FASTA_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/fasta_2/ # Path to the directory containing the input FASTA files.
OUTPUT_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/ # Please replace YourPath with the exact path on your machine # Path to the directory where the results will be saved.

for i in `ls $FASTA_DIR`; do
  echo $FASTA_DIR$i
   bash /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/alphafold-2.3.1/run_alphafold.sh \
  --fasta_paths=$FASTA_DIR$i \
  --max_template_date=2024-12-02 \
  --data_dir=$DATA_DIR \
  --output_dir=$OUTPUT_DIR \
  --model_preset=monomer_ptm \
#   --docker_image_name=alphafold_v3 \
  --models_to_relax=best \
  --db_preset=reduced_dbs \
  --num_multimer_predictions_per_model=1

#### Running an test alphafold in slurm
#!/bin/bash
#SBATCH --job-name=gpu_test_job    # Job name
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --time=02:00:00            # Max runtime: 2 hrs (adjust if needed)
#SBATCH --mem=50G                   # Memory request
#SBATCH --output=gpu_step2_output.log  # Output log file
#SBATCH --error=gpu_step2_error.log  # Error log file

bash step2_2.sh

