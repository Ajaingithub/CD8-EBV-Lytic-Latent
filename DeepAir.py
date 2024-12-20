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

# it is better to download the ProtTrans rather than using Pytorch
import torch
print(torch.__version__)  
print(torch.cuda.is_available())  
print(torch.cuda.get_device_name(0)) 

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

## Using GPUs
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

### Since we have to make some changes in the model.py. We cannot use the singularity image.
# I have used Kaliana lab to perform local alpha fold installation
# https://github.com/kalininalab/alphafold_non_docker
# cd /diazlab/data3/.abhinav/tools/miniconda3/envs/deepair/lib/python3.8/site-packages/ && patch -p0 < $alphafold_path/docker/openmm.patch

### Please keep all the python version compatible other you will fall into the python dependency hell.
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
# Available devices: [gpu(id=0)]

# Create a matrix of ones and perform a simple operation
x = jnp.ones((1000, 1000))  # Create a 1000x1000 matrix of ones
y = jnp.dot(x, x)           # Matrix multiplication (dot product)

# Print a small part of the result to confirm
print("Result of matrix multiplication (partial):", y[:5, :5])  

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
source activate alphafold2           # Use your specific Conda environment

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

# bash run_alphafold.sh -d /diazlab/data3/abhinav/resource/alphafold_database/ \
# -o /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/new_af/alphafold-2.3.1/output_af2_structure/ \
# -f /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/fasta_2/QQYDSHLYT_9680.fasta \
# -t 2020-05-16

# bash run_alphafold.sh -d /diazlab/data3/abhinav/resource/alphafold_database/ \
# -o /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/new_af/alphafold-2.3.1/output_af2_return_representation/ \
# -f /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/fasta_2/VYGSGSPSNWFHP_460.fasta \
# -t 2020-05-15

# bash run_alphafold.sh -d /diazlab/data3/abhinav/resource/alphafold_database/ \
# -o /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/new_af/alphafold-2.3.1/output_af2_return_representation/ \
# -f /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/fasta_2/QQYDNWPPWT_460.fasta \
# -t 2020-05-16

# End of script
# Running Step2
DATA_DIR=/diazlab/data3/abhinav/resource/alphafold_database/   # Path to the directory containing the AlphaFold 2 downloaded data.
FASTA_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/fasta_2/ # Path to the directory containing the input FASTA files.
OUTPUT_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/ 

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
#SBATCH --time=07-00:00:00            # Max runtime: 2 hrs (adjust if needed)
#SBATCH --mem=100G                   # Memory request
#SBATCH --output=gpu_step2_output.log  # Output log file
#SBATCH --error=gpu_step2_error.log  # Error log file

bash step2_2.sh

#### After completion of Step 2, we got the feature.pkl as well as the best model performance from af2
## I have made some changes in the step3 since it is picking up the relaxed_model_ptm. In af2 I could not fina any model with ptm instead there is monomer ptm
# so changed the name from result_model_1_ptm_pred_0.pkl to result_model_1_pred_0.pkl

# python step3.py \
# --AIR_file_path /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/QQYDNWPPWT_460_CDR3.csv \
# --AF2_feature_folder /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/new_af/alphafold-2.3.1/output_af2_return_representation/ \
# --output_folder /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/QQYDNWPPWT_460_output_feature

### After generating the files lets run DeepAIR
(2) For binding affinity prediciton (BAP) (Regression)
python ./maincode/DeepAIR_BAP.py  \
--input_data_file /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/sampledata/BAP/A0301_KLGGALQAK_IE-1_CMV_Reg.csv  \
--result_folder /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/sampledata/DeepAir_result \
--epitope A0301_KLGGALQAK_IE-1_CMV  \
--AF2_feature_folder /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/sampledata/structure_feature/BAP/10X  \
--transformer_model_folder /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/ProtTrans/prot_bert_bfd



##### After Performing the test run on the tools data. Testing it on our EBV CD8 test data at this location
# /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data
# while preparing the AIR file path we have kept the same header, as a default parameter 

Step1
conda activate DeepAIR
python /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/step1.py \
    --AIR_file_path /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/CD8_Ag_test_CDR3.csv \
    --output_table /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/CD8_Ag_test_CDR3Region.csv \
    --output_fasta_folder /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/fasta/

Step2: Running Alphafold
conda activate alphafold2
cd /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/new_af/alphafold-2.3.1/
bash run_alphafold.sh \
  -d /diazlab/data3/abhinav/resource/alphafold_database/ \
    -o /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/af2_output_structure/ \
      -f  /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/fasta/CVVNMGTLTF_1.fasta \
        -t 2020-05-16


### Running on the GPU 
#!/bin/bash
#SBATCH --job-name=gpu_test_job    # Job name
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --time=02:00:00            # Max runtime: 2 hrs
#SBATCH --mem=8G                   # Memory request
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

# End of script
# Running Step2
source ~/.bashrc
source /diazlab/data3/.abhinav/tools/miniconda3/etc/profile.d/conda.sh
source activate alphafold2
conda init
conda activate alphafold2

DATA_DIR=/diazlab/data3/abhinav/resource/alphafold_database/   # Path to the directory containing the AlphaFold 2 downloaded data.
FASTA_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/fasta/fasta2/ # Path to the directory containing the input FASTA files.
OUTPUT_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/af2_output_structure/ 

cd /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/new_af/alphafold-2.3.1/
for i in `ls $FASTA_DIR`; do
  echo $FASTA_DIR$i
   bash /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/alphafold-2.3.1/run_alphafold.sh \
   -d /diazlab/data3/abhinav/resource/alphafold_database/ \
    -o /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/af2_output_structure/ \
    -f $FASTA_DIR$i \
    -t 2023-05-16

#### Running an test alphafold in slurm
#!/bin/bash
#SBATCH --job-name=gpu_test_job    # Job name
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --time=07-00:00:00            # Max runtime: 2 hrs (adjust if needed)
#SBATCH --mem=100G                   # Memory request
#SBATCH --output=gpu_step2_output.log  # Output log file
#SBATCH --error=gpu_step2_error.log  # Error log file

bash step2_2.sh


#!/bin/bash
#SBATCH --job-name=test_multiple_job    # Job name
#SBATCH --gres=gpu:1                    # Request 1 GPU
#SBATCH --time=02:00:00                 # Max runtime: 2 hrs (adjust if needed)
#SBATCH --mem=32G                        # Memory request
#SBATCH --output=test_multiple_output_%A_%a.log  # Output log file
#SBATCH --error=gpu_step2_error_%A_%a.log  # Error log file
#SBATCH --array=0-$(($(ls $FASTA_DIR | wc -l) - 1))  # Array of jobs, one for each FASTA file

source ~/.bashrc
source /diazlab/data3/.abhinav/tools/miniconda3/etc/profile.d/conda.sh
source activate alphafold2
conda init
conda activate alphafold2

# Directories
DATA_DIR=/diazlab/data3/abhinav/resource/alphafold_database/
FASTA_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/fasta/fasta2/
OUTPUT_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/data/af2_output_structure/

# Get the list of FASTA files
FASTA_FILES=($(ls $FASTA_DIR))

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

