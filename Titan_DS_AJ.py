import argparse
import json
import logging
import os
import sys

import numpy as np
import pandas as pd
from pytoda.files import read_smi
from sklearn.metrics import (
    accuracy_score, auc, average_precision_score, balanced_accuracy_score,
    confusion_matrix, matthews_corrcoef, precision_score, roc_curve
)

from paccmann_tcr.models import knn
from paccmann_tcr.utils.plot_knn import plot_roc_prc
from paccmann_tcr.utils.utils import cutoff_youdens_j

train_epi = read_smi("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/epitopes.smi", names=['data'])
train_tcr = read_smi("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/tcr_full.csv", names=['data'])

test_epitope_fp = "."
test_tcr_fp = "."
test_epi = train_epi if test_epitope_fp == '.' else read_smi(
    test_epitope_fp, names=['data']
    )
test_tcr = train_tcr if test_tcr_fp == '.' else read_smi(
    test_tcr_fp, names=['data']
    )

import os
import pandas as pd

# Assuming `data_path`, `train_name`, `test_name`, `train_epi`, `train_tcr`, `test_epi`, `test_tcr`, and `logger` are already defined
data_path = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/strict_split/"
train_name = "train+covid.csv"
test_name  = "test+covid.csv"

# Load train and test data for the given fold
fold = 0
train_data = pd.read_csv(
    os.path.join(data_path, f'fold{fold}', train_name), index_col=0
)
test_data = pd.read_csv(
    os.path.join(data_path, f'fold{fold}', test_name), index_col=0
)

# Extract epitopes and TCRs for the training and testing datasets
train_epitopes = train_epi.loc[train_data['ligand_name']]['data']
train_tcrs = train_tcr.loc[train_data['sequence_id']]['data']
test_epitopes = test_epi.loc[test_data['ligand_name']]['data']
test_tcrs = test_tcr.loc[test_data['sequence_id']]['data']

# Extract labels for training and testing data
train_labels = train_data['label']
test_labels = test_data['label']


### Running KNN
# predictions, knn_labels = knn(
#     train_epitopes,
#     train_tcrs,
#     train_labels,
#     test_epitopes,
#     test_tcrs,
#     k=25,
#     return_knn_labels=True
# )

from typing import Iterable, List
import numpy as np
from Levenshtein import distance as levenshtein_distance
k = 25
test_epitopes.iloc[0]
test_tcrs.iloc[0]
predictions, knn_labels = [], []
for epitope, tcr in zip(test_epitopes, test_tcrs):
    # test_epitopes.iloc[0]
    # test_tcrs.iloc[0]
    el = len(epitope)
    tl = len(tcr)
    # Calculate Levenshtein distances for epitopes and TCRs, normalized by their lengths
    epitope_dists = [levenshtein_distance(epitope, e) / el for e in train_epitopes]
    tcr_dists = [levenshtein_distance(tcr, t) / tl for t in train_tcrs]
    # Find indices of k-nearest neighbors based on combined distances
    knns = np.argsort(np.array(epitope_dists) + np.array(tcr_dists))[:k]
    # Get the labels of the k-nearest neighbors
    _knn_labels = np.array(train_labels)[knns]
    # Append the mean label of k-nearest neighbors to predictions
    predictions.append(np.mean(_knn_labels))
    # Optionally, store the actual labels of the k-nearest neighbors
    knn_labels.append(_knn_labels)


# python3 knn_cv.py \
# -d /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/tcr_split \
# -tr train+covid.csv \
# -te test+covid.csv \
# -f 10 \
# -ep /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/epitopes.csv \
# -tcr /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/tcr_full.csv \
# -r /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/test_result/ \
# -k 25

# Transformers Attention 
## Directly running
# python3 ../scripts/flexible_training.py \
# data/train_small.csv \
# data/test_small.csv  \
# data/tcr_full.csv \
# data/epitopes.smi \
# trained_model \
# data/params_small.json \
# tutorial_setting \
# bimodal_mca

# Import packages
import os
import numpy as np
import torch
import json

## it is pattern utility to create machine learning models. It is to dynamically initialize and retrieve models supported by PaccMann predictor framework.
# model_config = {
#     'model_name': 'mlp',  # Specify the model type (e.g., MLP, RNN, etc.)
#     'input_dim': 512,     # Input dimension
#     'hidden_dim': 256,    # Hidden dimension
#     'output_dim': 1       # Output dimension (e.g., binary classification)
# }
# model = MODEL_FACTORY[model_config['model_name']](**model_config)
from paccmann_predictor.models import MODEL_FACTORY 
from paccmann_predictor.utils.utils import get_device ### cuda if GPU, otherwise CPU

## DrugAffinityDataset: A dataset class for handling drug-target interaction data. Each instance typically contains:
#     A drug molecule (in SMILES format).
#     A target protein (sequence or embedding).
#     A binding affinity or interaction score.
# ProteinProteinInteractionDataset: A dataset class for protein-protein interaction data. Each instance typically contains:
#     Two protein sequences.
#     A label indicating interaction (e.g., 0 for no interaction, 1 for interaction).
# Why they’re useful: They handle data loading, tokenization, and preprocessing out of the box, allowing you to focus on model building and training.

from pytoda.datasets import (
    DrugAffinityDataset, ProteinProteinInteractionDataset
)

# 4. ProteinFeatureLanguage and ProteinLanguage from pytoda.proteins
#         ProteinLanguage: A class for handling protein sequences as language models. It tokenizes protein sequences into smaller units (amino acids or sub-sequences) 
# and represents them numerically for input into models.
#         ProteinFeatureLanguage: A more advanced version that includes additional features (e.g., secondary structure, chemical properties).
# Why they’re useful: Proteins are complex sequences, and these classes handle tokenization, padding, and feature extraction, which are necessary for deep learning models.
from pytoda.proteins import ProteinFeatureLanguage, ProteinLanguage

# What it is: A tokenizer specifically for SMILES (Simplified Molecular Input Line Entry System), a representation of chemical structures.
# Purpose: Converts SMILES strings into tokenized forms that machine learning models can understand.
# Use case:
from pytoda.smiles.smiles_language import SMILESTokenizer
# Tokenization is basically breaking down into k-mer sequences like
# Input: "MVLSPADKTNVKAA"
# Tokens: ['MVL', 'VLS', 'LSP', 'SPA', 'PAD', 'ADK', 'DKT', 'KTN', 'TNV', 'NVK', 'VKA', 'KAA']

from sklearn.metrics import (
    auc, average_precision_score, precision_recall_curve, roc_curve
)

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm

device = get_device()

### Loading the data
data_dir = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/tutorial/data/"
tcrs = os.path.join(data_dir + "tcr_full.csv")
epitopes = os.path.join(data_dir, 'epitopes.smi')
test_data = os.path.join(data_dir, 'test_small.csv')
# model_path = os.path.join('/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/trained_model/')
model_path = os.path.join('/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/tutorial/data/')

# Process parameter file:
# params_filepath = os.path.join(model_path, 'model_params.json')
params_filepath = os.path.join(model_path, 'params_small.json')
params = {}
with open(params_filepath) as fp:
    params.update(json.load(fp))

# Load languages
# The tokenizer contains predefined rules for breaking down SMILES strings into tokens (e.g., [C], [O], or [=]) for input to a neural network.
# at this location /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/trained_model/token_count.json
smiles_language = SMILESTokenizer.from_pretrained(model_path)

# The set_encoding_transforms method configures how SMILES strings will be processed before being passed to the model:
#     randomize=None:
#         No randomization of SMILES tokens (e.g., skipping augmentation that generates multiple valid SMILES for the same molecule).
#     add_start_and_stop=params.get('ligand_start_stop_token', True):
#         Adds special "start" and "stop" tokens to the beginning and end of each SMILES string, as specified in the params dictionary (True by default).
#         Useful for models using sequence-based architectures (e.g., RNNs, Transformers).
#     padding=params.get('ligand_padding', True):
#         Enables padding to ensure that all SMILES sequences are of the same length in a batch.
#     padding_length=params.get('ligand_padding_length', True):
#         Sets the maximum length for padding sequences, defined in the model parameters (500 in this case).

smiles_language.set_encoding_transforms(
    randomize=None,
    add_start_and_stop=params.get('ligand_start_stop_token', True),
    padding=params.get('ligand_padding', True),
    padding_length=params.get('ligand_padding_length', True),
)

smiles_language.set_smiles_transforms(
    augment=False,
)

# ProteinFeatureLanguage.load:
#     Loads a pre-trained protein feature representation (e.g., BLOSUM matrix) from a file (protein_language.pkl).
#     This is used for predefined embeddings.

if params.get('receptor_embedding', 'learned') == 'predefined':
    protein_language = ProteinFeatureLanguage.load(
        os.path.join(model_path, 'protein_language.pkl')
    )
else:
    protein_language = ProteinLanguage.load(
        os.path.join(model_path, 'protein_language.pkl')
    )

# Load the data
# test_dataset = DrugAffinityDataset(
#     drug_affinity_filepath=test_data,
#     smi_filepath=epitopes,
#     protein_filepath=tcrs,
#     smiles_language=smiles_language,
#     protein_language=protein_language,
#     smiles_padding=params.get('ligand_padding', True),
#     smiles_padding_length=params.get('ligand_padding_length', None),
#     smiles_add_start_and_stop=params.get(
#         'ligand_add_start_stop', True
#     ),
#     protein_padding=params.get('receptor_padding', True),
#     protein_padding_length=params.get('receptor_padding_length', None),
#     protein_add_start_and_stop=params.get(
#         'receptor_add_start_stop', True
#     ),
#     drug_affinity_dtype=torch.float,
#     backend='eager',
#     iterate_dataset=True
# )

test_dataset = DrugAffinityDataset(
    drug_affinity_filepath=test_data,
    smi_filepath=epitopes,
    protein_filepath=tcrs,
    smiles_language=smiles_language,
    protein_language=protein_language,
    drug_affinity_dtype=torch.float,
    backend='eager',
    iterate_dataset=True
)


test_loader = torch.utils.data.DataLoader(
    dataset=test_dataset,
    batch_size=10,
    shuffle=False,
    drop_last=False,
    num_workers=params.get('num_workers', 0)
)
