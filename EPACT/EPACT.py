### After installing EPACT.py
### I have made a small change in the source file plm_trainer.py
### in line 20
# class EpitopeLMTrainer(BaseTrainer):
#     def __init__(self, config):
#         super(EpitopeLMTrainer, self).__init__(config)
        
#         # Check if CUDA is available
#         if torch.cuda.is_available():
#             # You can choose a specific device like 'cuda:0' or 'cuda:1'
#             self.device = torch.device("cuda:0")  # Change 0 to the desired GPU index
#         else:
#             self.device = torch.device("cpu")  # Use CPU if CUDA is not available

#         # Initialize model and move to the selected device
#         self.model = PeptideLM(**config.model)
#         self.model.to(self.device)  # Move model to the device (CUDA or CPU)

# conda activate EPACT_env2

from argparse import ArgumentParser
import os

import numpy as np
import pandas as pd
from torch.utils.data import DataLoader, random_split

os.chdir("/diazlab/data3/.abhinav/tools/")
from EPACT.utils import set_seed, load_config
from EPACT.dataset import EpitopeDataset, PairedCDR3Dataset, PairedCDR123Dataset, EpitopeBatchConverter, PairedCDR3BatchConverter, PairedCDR123BatchConverter
from EPACT.trainer import EpitopeLMTrainer, PairedTCRLMTrainer
os.chdir("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/")
os.makedirs("EPACT", exist_ok=True)
os.chdir("./EPACT/")

tool_path = "/diazlab/data3/.abhinav/tools/EPACT"
config_path = os.path.join(tool_path,"configs/config-pretrain-epitope-lm.yml")
config = load_config(config_path)
    
if not os.path.exists(config.training.log_dir):
    os.makedirs(config.training.log_dir)

os.system(f'cp {config_path} {config.training.log_dir}/config.yml')

set_seed(config.training.seed)
    
import torch
print(torch.cuda.device_count())  # Number of available GPUs
print(torch.cuda.current_device())  # Current device index
print(torch.cuda.get_device_name(0))  # Name of GPU at index 0 (if available)
# device = torch.device("cuda:0")  # Ensure 0 is a valid device index

print(torch.version.cuda)  # Check CUDA version used by PyTorch
print(torch.backends.cudnn.version())  # Check cuDNN version

if config.task == 'epitope-lm':
    dataset = EpitopeDataset(os.path.join(tool_path,config.data.data_path))
    train_size = int(0.8 * len(dataset))  # 80% for training
    val_size = int(0.1 * len(dataset))    # 10% for validation
    test_size = len(dataset) - train_size - val_size  # Remaining for testing
    # Perform random split
    train_data, val_data, test_data = random_split(dataset, [train_size, val_size, test_size])
    batch_converter = EpitopeBatchConverter(max_epitope_len=config.data.max_epi_len, use_atchley_factor=True)
    Trainer = EpitopeLMTrainer(config)
    
elif config.task == 'cdr3-lm':
        df = pd.read_csv(config.data.data_path)
        dataset = PairedCDR3Dataset(config.data.data_path)
        train_size = int(0.8 * len(dataset))  # 80% for training
        val_size = int(0.1 * len(dataset))    # 10% for validation
        test_size = len(dataset) - train_size - val_size  # Remaining for testing
        # Perform random split
        train_data, val_data, test_data = random_split(dataset, [train_size, val_size, test_size])
        batch_converter = PairedCDR3BatchConverter(max_cdr3_len=config.data.max_cdr3_len, use_atchley_factor=True, mask_prob=0.15, beta_mask_prob=0.2)
        Trainer = PairedTCRLMTrainer(config)
    
elif config.task == 'cdr123-lm':
    dataset = PairedCDR123Dataset(config.data.data_path)    
    train_size = int(0.8 * len(dataset))  # 80% for training
    val_size = int(0.1 * len(dataset))    # 10% for validation
    test_size = len(dataset) - train_size - val_size  # Remaining for testing
    # Perform random split
    train_data, val_data, test_data = random_split(dataset, [train_size, val_size, test_size])        
    batch_converter = PairedCDR123BatchConverter(max_cdr3_len=config.data.max_cdr3_len, use_atchley_factor=True, mask_prob=0.1, beta_mask_prob=0.2)
    Trainer = PairedTCRLMTrainer(config)

print(f'Train data size: {len(train_data)}')
print(f'Val data size: {len(val_data)}')
print(f'Test data size: {len(test_data)}')
    
train_loader = DataLoader(train_data, batch_size=config.training.batch_size, collate_fn=batch_converter,num_workers=config.training.num_workers, shuffle=True)
val_loader = DataLoader(val_data, batch_size=config.training.batch_size, collate_fn=batch_converter,num_workers=config.training.num_workers, shuffle=False)
test_loader = DataLoader(test_data, batch_size=config.training.batch_size, collate_fn=batch_converter, num_workers=config.training.num_workers, shuffle=False)

Trainer.fit(train_loader, val_loader)
Trainer.test(test_loader, model_location=os.path.join(config.training.log_dir, 'checkpoint.pt'))
