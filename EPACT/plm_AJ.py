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
