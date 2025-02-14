#!/bin/bash
#SBATCH --job-name=epitope_train    # Job name
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --time=10:00:00            # Max runtime: 10 hrs
#SBATCH --mem=16G                   # Memory request
#SBATCH --output=epitope_train.log  # Output log file
#SBATCH --error=epitope_train.log  # Error log file

# Activate your Conda environment
source ~/.bashrc
source /diazlab/data3/.abhinav/tools/miniconda3/etc/profile.d/conda.sh
source activate EPACT_env2
conda init
conda activate EPACT_env2

cd /diazlab/data3/.abhinav/tools/EPACT/
# pretrain epitope masked language model.
python scripts/pretrain/pretrain_plm.py --config configs/config-pretrain-epitope-lm.yml

# pretrain paired cdr3 masked language model.
# python scripts/pretrain/pretrain_plm.py --config configs/config-pretrain-cdr3-lm.yml

# pretrain paired cdr123 masked language model.
# python scripts/pretrain/pretrain_plm.py --config configs/config-pretrain-cdr123-lm.yml