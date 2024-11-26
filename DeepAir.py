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
