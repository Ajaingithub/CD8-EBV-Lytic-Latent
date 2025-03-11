import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import random

# Step 1: Define the Tokenizer
class ProteinTokenizer:
    def __init__(self):
        self.amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        self.token2idx = {aa: idx + 1 for idx, aa in enumerate(self.amino_acids)}
        self.token2idx["PAD"] = 0  # Padding token
        self.token2idx["MASK"] = len(self.token2idx)  # Mask token
        self.idx2token = {idx: aa for aa, idx in self.token2idx.items()}
    
    def encode(self, sequence):
        return [self.token2idx[aa] for aa in sequence]
    
    def decode(self, tokens):
        return "".join([self.idx2token[idx] for idx in tokens if idx in self.idx2token])

tokenizer = ProteinTokenizer()

# Step 2: Define the Transformer Model
class TCRLM(nn.Module):
    def __init__(self, num_layers=3, embed_dim=128, num_heads=4, max_seq_length=20, dropout=0.1):
        super(TCRLM, self).__init__()
        self.max_seq_length = max_seq_length
        self.embed_dim = embed_dim
        self.token_embedding = nn.Embedding(len(tokenizer.token2idx), embed_dim)
        self.position_embedding = nn.Embedding(max_seq_length, embed_dim)
        self.transformer = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(embed_dim, num_heads, dim_feedforward=256, dropout=dropout),
            num_layers=num_layers
        )
        self.lm_head = nn.Linear(embed_dim, len(tokenizer.token2idx))
        self.dropout = nn.Dropout(dropout)
    
    def forward(self, x):
        seq_length = x.size(1)
        position_ids = torch.arange(seq_length, dtype=torch.long, device=x.device).unsqueeze(0)
        x = self.token_embedding(x) + self.position_embedding(position_ids)
        x = self.transformer(x)
        x = self.lm_head(self.dropout(x))
        return x

# Step 3: Create a Dataset for Training
class ProteinDataset(Dataset):
    def __init__(self, sequences):
        self.sequences = sequences
        self.max_seq_length = 20
        self.mask_token = tokenizer.token2idx["MASK"]
    
    def __len__(self):
        return len(self.sequences)
    
    def __getitem__(self, idx):
        seq = self.sequences[idx]
        tokenized = tokenizer.encode(seq)
        if len(tokenized) < self.max_seq_length:
            tokenized += [0] * (self.max_seq_length - len(tokenized))
        masked_seq = tokenized.copy()
        mask_idx = random.randint(0, len(seq) - 1)
        masked_seq[mask_idx] = self.mask_token  # Apply masking
        return torch.tensor(masked_seq), torch.tensor(tokenized)

# Sample TCR sequences
tcr_sequences = ["CASSL", "CASSP", "CAWSV", "CSARD", "CATSD"]
dataset = ProteinDataset(tcr_sequences)
dataloader = DataLoader(dataset, batch_size=2, shuffle=True)

# Step 4: Train the Model
model = TCRLM()
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

def train(model, dataloader, epochs=10):
    model.train()
    for epoch in range(epochs):
        total_loss = 0
        for masked_seq, target_seq in dataloader:
            optimizer.zero_grad()
            output = model(masked_seq)
            loss = criterion(output.view(-1, len(tokenizer.token2idx)), target_seq.view(-1))
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        print(f"Epoch {epoch + 1}, Loss: {total_loss:.4f}")

train(model, dataloader)

# while on the c4-n38 and c4-n39 it is NVIDIA A100 80GB PCIe MIG 1g.10gb i.e. it is divided into 7 10GB MIGs. While the head node GPU node NVIDIA RTX A4000 use the full
# memory to run which is much faster. We cannot use 2 instances at a time which makes it very slow. So the admin has to provide the entire GPU to train for a time being. 

#### Predicting the TCR-pMHC using the EPACT dataset
# predict cross-validation results
cd /diazlab/data3/.abhinav/tools/EPACT
for i in {1..5}
do
    python scripts/predict/predict_tcr_pmhc_binding.py \
        --config configs/config-paired-cdr123-pmhc-binding.yml \
        --input_data_path data/binding/Full-TCR/k-fold-data/val_fold_${i}.csv \
        --model_location checkpoints/paired-cdr123-pmhc-binding/paired-cdr123-pmhc-binding-model-fold-${i}.pt\
        --log_dir results/preds-cdr123-pmhc-binding/Fold_${i}/
done

# predict distance matrices and contact sites between MEL8 TCR and HLA-A2-presented peptides.
cd /diazlab/data3/.abhinav/tools/EPACT
for i in {1..5}
do
    python scripts/predict/predict_tcr_pmhc_interact.py --config configs/config-paired-cdr123-pmhc-interact.yml \
        --input_data_path data/MEL8_A0201_peptides.csv \
        --model_location checkpoints/paired-cdr123-pmhc-interaction/paired-cdr123-pmhc-interaction-model-fold-${i}.pt \
        --log_dir results/interaction-MEL8-bg-cdr123-closest/Fold_${i}/
done

import pickle

# Open the .pkl file in read-binary mode
with open('/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/EPACT/results/interaction-CD8-EBV-bg-cdr123-closest_all_unique/predictions.pkl', 'rb') as file:
    # Load the object from the file
    data = pickle.load(file)

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import pandas as pd
savedir = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/EPACT/results/interaction-CD8-EBV-bg-cdr123-closest_all_unique/"

# Collect the distance matrices for epitope "ALWALPHAA"
epitopes = list(set([entry['epitope'] for entry in data]))
cdrs = ["cdr1.alpha","cdr3.alpha","cdr3.beta"]

for epit in epitopes:
    for cdr in cdrs:
        distance_matrices = []
        for entry in data:
            if entry['epitope'] == epit:
                distance_matrices.append(entry['dist'][cdr])
        
        concatenated_matrix = np.concatenate(distance_matrices, axis=0)
        amino_acids = list(epit)
        column_means = np.mean(concatenated_matrix, axis=0)
        column_std = np.std(concatenated_matrix, axis=0)

        plt.figure(figsize=(5, 3))
        bars = plt.bar(range(len(column_means)), column_means, yerr=column_std, capsize=5, color='skyblue', label='Mean Values')
        plt.xticks(range(len(amino_acids)), amino_acids)
        plt.xlabel("Epitope "+epit+" amino acid")
        plt.ylabel(cdr +" distance")
        plt.title(epit)
        plt.legend()
        plt.tight_layout()
        os.makedirs(Path(savedir,"distance/"), exist_ok=True)
        full_path = Path(savedir,"distance/") / (epit + "_" + cdr +"_distance.pdf")
        plt.savefig(full_path)
        df = pd.DataFrame(data=concatenated_matrix, columns = amino_acids)  # 1st row as the column names
        full_path = Path(savedir,"distance/") / (epit + "_" + cdr +"_distance.csv")
        df.to_csv(full_path, index =  False)


#### Contact
epitopes = list(set([entry['epitope'] for entry in data]))
cdrs = ["cdr1.alpha","cdr3.alpha","cdr3.beta"]

for epit in epitopes:
    for cdr in cdrs:
        distance_matrices = []
        for entry in data:
            if entry['epitope'] == epit:
                distance_matrices.append(entry['contact'][cdr])
        
        concatenated_matrix = np.concatenate(distance_matrices, axis=0)
        amino_acids = list(epit)
        column_means = np.mean(concatenated_matrix, axis=0)
        column_std = np.std(concatenated_matrix, axis=0)

        plt.figure(figsize=(5, 3))
        bars = plt.bar(range(len(column_means)), column_means, yerr=column_std, capsize=5, color='forestgreen', label='Mean Values')
        plt.xticks(range(len(amino_acids)), amino_acids)
        plt.xlabel("Epitope "+epit+" amino acid")
        plt.ylabel(cdr +" contact score")
        plt.title(epit)
        plt.legend()
        plt.tight_layout()
        os.makedirs(Path(savedir,"contact/"), exist_ok=True)
        full_path = Path(savedir,"contact/") / (epit + "_" + cdr +"_contact.pdf")
        plt.savefig(full_path)
        df = pd.DataFrame(data=concatenated_matrix, columns = amino_acids)  # 1st row as the column names
        full_path = Path(savedir,"contact/") / (epit + "_" + cdr +"_contact.csv")
        df.to_csv(full_path, index =  False)

