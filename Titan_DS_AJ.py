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
print(device)

### Loading the data
data_dir = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/tutorial/data/"
tcrs = os.path.join(data_dir + "tcr_full.csv")
epitopes = os.path.join(data_dir, 'epitopes.smi')
test_data = os.path.join(data_dir, 'test_small.csv')
# model_path = os.path.join('/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/trained_model/')
model_path = os.path.join('/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/tutorial/trained_model/tutorial_setting/')
# model_path = os.path.join('/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/tutorial/TITAN_small/')

# This setup configures a convolutional + attention-based model for ligand-receptor interaction prediction. It processes ligands and receptors via embedding layers, 
# applies convolutional layers for feature extraction, focuses on important regions via attention, and combines the extracted features using dense layers. The 
# regularization (dropout) and batch size are tuned for efficient training.
### Model Parameters
{"augment_smiles": false, "smiles_canonical": false, "ligand_start_stop_token": true, "receptor_start_stop_token": true, "ligand_padding_length": 422, 
 "receptor_padding_length": 130, "receptor_embedding": "learned", "ligand_embedding": "learned", "predefined_embedding": "blosum", "receptor_embedding_size": 26, 
 "ligand_embedding_size": 26, "ligand_filters": [26, 26], "receptor_filters": [26, 26], "ligand_kernel_sizes": [[5, 26], [7, 26]], "receptor_kernel_sizes": [[3, 26], 
 [5, 26]], "ligand_attention_size": 32, "receptor_attention_size": 32, "dense_hidden_sizes": [200, 80], "activation_fn": "relu", "dropout": 0.5, "batch_norm": false, 
 "batch_size": 512, "lr": 0.0001, "epochs": 100, "save_model": 20, "ligand_vocabulary_size": 575, "receptor_vocabulary_size": 33, "ligand_as": "smiles", 
 "number_of_parameters": 418369}
# receptor_embedding: "learned"`
# The receptor's input (e.g., a protein sequence) is converted into embeddings that are learned during training. This means the model will optimize these embeddings 
# as part of its parameters.
# ligand_embedding: "learned"`
# Similarly, the ligand's input (e.g., a SMILES string) is also converted into embeddings that are learned during training.
# predefined_embedding: "blosum"`
# This specifies a predefined embedding (e.g., BLOSUM matrix for amino acids). The model may use it as a static embedding for receptors or as initialization for 
# learned embeddings.

# 2. Embedding Sizes
# These define the dimensions of the embeddings (vectors representing ligands and receptors):
#     receptor_embedding_size: 26
#     Each receptor token (e.g., an amino acid) is represented as a 26-dimensional vector.
#     Example: A -> [0.1, 0.3, ..., -0.2] (26 values).
#     ligand_embedding_size: 26
#     Each ligand token (e.g., a SMILES character) is represented similarly as a 26-dimensional vector.

# 3. Convolutional Filters
# These define the number of filters (feature maps) and filter sizes (kernel sizes) used in convolutional layers for feature extraction from ligands and receptors.
#     ligand_filters: [26, 26]
#     The number of convolutional filters applied to ligand embeddings at each layer. Here, two layers each with 26 filters.
#     receptor_filters: [26, 26]
#     Similarly, the number of filters for receptor embeddings.
#     ligand_kernel_sizes: [[5, 26], [7, 26]]
#     The size of the convolutional kernels for ligand sequences.
#         The first layer uses a 5x26 kernel (operates over a window of 5 tokens).
#         The second layer uses a 7x26 kernel.
#     receptor_kernel_sizes: [[3, 26], [5, 26]]
#     Kernel sizes for receptor sequences.
#         The first layer uses a 3x26 kernel.
#         The second layer uses a 5x26 kernel.

# 4. Attention Parameters
# These control the attention mechanisms, which allow the model to focus on the most relevant parts of the ligand or receptor sequence:
#     ligand_attention_size: 32
#     The size of the attention layer for ligands. This defines the attention weight vector's dimension.
#     receptor_attention_size: 32
#     The size of the attention layer for receptors.

# 5. Dense Layers
# These parameters define the dense (fully connected) layers that combine features extracted from ligands and receptors:
#     dense_hidden_sizes: [200, 80]
#     The dense layer structure:
#         First dense layer: 200 neurons.
#         Second dense layer: 80 neurons.
#     **activation_fn: "relu"** The activation function used in the dense layers. ReLU` (Rectified Linear Unit) introduces non-linearity.

# 6. Regularization Parameters
# These prevent overfitting and stabilize training:
#     dropout: 0.5
#     Dropout randomly disables 50% of neurons during training to prevent overfitting.
#     batch_norm: false
#     Whether to apply batch normalization (standardizing intermediate activations). Here, it’s disabled.

# 7. Training Parameters
# These define how the model processes data during training:
#     batch_size: 512
#     The number of samples processed together in one forward/backward pass during training.

# Process parameter file:
# params_filepath = os.path.join(model_path, 'model_params.json')
params_filepath = os.path.join(model_path, 'model_params.json')
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
test_dataset = DrugAffinityDataset(
    drug_affinity_filepath=test_data,
    smi_filepath=epitopes,
    protein_filepath=tcrs,
    smiles_language=smiles_language,
    protein_language=protein_language,
    smiles_padding=params.get('ligand_padding', True),
    smiles_padding_length=params.get('ligand_padding_length', None),
    smiles_add_start_and_stop=params.get(
        'ligand_add_start_stop', True
    ),
    protein_padding=params.get('receptor_padding', True),
    protein_padding_length=params.get('receptor_padding_length', None),
    protein_add_start_and_stop=params.get(
        'receptor_add_start_stop', True
    ),
    drug_affinity_dtype=torch.float,
    backend='eager',
    iterate_dataset=True
)

# dataset=test_dataset
# The test_dataset object is passed as the data source. This can be any PyTorch Dataset object (like ProteinSequenceDataset in your case). The DataLoader will pull batches of data from this dataset.
# batch_size=10
# Specifies the number of samples to include in each batch.
#     Each time you call the DataLoader, it will return a batch of 10 samples from the dataset.

# shuffle=False
# Controls whether the data is shuffled before being loaded.
#     False: The data is read sequentially (in the order defined in the dataset).
#     True: The data order is randomized, which is common for training but not needed during testing or inference.

# drop_last=False
# Determines whether to drop the last incomplete batch if the dataset size is not divisible by the batch size.
#     False: Include the smaller, last batch (if the dataset size isn’t perfectly divisible by batch_size).
#     True: Discard the smaller last batch.

# num_workers=params.get('num_workers', 0)
# Specifies the number of subprocesses used for data loading.
#     0: Data loading is performed in the main process (slower, but no multiprocessing overhead).
#     >0: Enables parallel data loading using multiple worker processes, speeding up loading for large datasets.
#     The value is retrieved from the params dictionary using .get(). If num_workers is not found in params, it defaults to 0.

test_loader = torch.utils.data.DataLoader(
    dataset=test_dataset,
    batch_size=10,
    shuffle=False,
    drop_last=False,
    num_workers=params.get('num_workers', 0)
)

# Load the Model
model_fn = params.get('model_fn', 'bimodal_mca')
model = MODEL_FACTORY[model_fn](params).to(device)
model
# BimodalMCA(
#   (loss_fn): BCELoss()
#   (act_fn): ReLU()
#   (ligand_embedding): Embedding(575, 26)
#   (receptor_embedding): Embedding(33, 26)
#   (ligand_convolutional_layers): Sequential(
#     (ligand_convolutional_0): Sequential(
#       (convolve): Conv2d(1, 26, kernel_size=(5, 26), stride=(1, 1), padding=(2, 0))
#       (squeeze): Squeeze()
#       (act_fn): ReLU()
#       (dropout): Dropout(p=0.5, inplace=False)
#       (batch_norm): Identity()
#     )
#     (ligand_convolutional_1): Sequential(
#       (convolve): Conv2d(1, 26, kernel_size=(7, 26), stride=(1, 1), padding=(3, 0))
#       (squeeze): Squeeze()
#       (act_fn): ReLU()
#       (dropout): Dropout(p=0.5, inplace=False)
#   (context_attention_receptor_layers): Sequential(
#     (context_attention_receptor_0): ContextAttentionLayer(
#       (individual_nonlinearity): Sequential()
#       (reference_projection): Sequential(
#         (projection): Linear(in_features=26, out_features=32, bias=True)
#         (act_fn): Sequential()
#       )
#       (context_projection): Sequential(
#         (projection): Linear(in_features=26, out_features=32, bias=True)
#         (act_fn): Sequential()
#       )
#       (context_hidden_projection): Sequential(
#         (projection): Linear(in_features=422, out_features=130, bias=True)
#         (act_fn): Sequential()
#       )
#       (alpha_projection): Sequential(
#         (projection): Linear(in_features=32, out_features=1, bias=False)
#         (squeeze): Squeeze()
#       (context_projection): Sequential(
#         (projection): Linear(in_features=26, out_features=32, bias=True)
#         (act_fn): Sequential()
#       )
#       (context_hidden_projection): Sequential(
#         (projection): Linear(in_features=422, out_features=130, bias=True)
#         (act_fn): Sequential()
#       )
#       (alpha_projection): Sequential(
#         (projection): Linear(in_features=32, out_features=1, bias=False)
#         (squeeze): Squeeze()
#         (temperature): Temperature()
#         (softmax): Softmax(dim=1)
#       )
#     )
#     (context_attention_receptor_2): ContextAttentionLayer(
#       (individual_nonlinearity): Sequential()
#       (reference_projection): Sequential(
#         (projection): Linear(in_features=26, out_features=32, bias=True)
#         (act_fn): Sequential()
#       )
#       (context_projection): Sequential(
#         (projection): Linear(in_features=26, out_features=32, bias=True)
#         (act_fn): Sequential()
#       )
#       (context_hidden_projection): Sequential(
#         (projection): Linear(in_features=422, out_features=130, bias=True)
#         (act_fn): Sequential()
#       )
#       (alpha_projection): Sequential(
#         (projection): Linear(in_features=32, out_features=1, bias=False)
#         (squeeze): Squeeze()
#         (temperature): Temperature()
#         (softmax): Softmax(dim=1)
#       )
#     )
#   )
#   (dense_layers): Sequential(
#     (dense_0): Sequential(
#       (projection): Linear(in_features=156, out_features=200, bias=True)
#       (batch_norm): Identity()
#       (act_fn): ReLU()
#       (dropout): Dropout(p=0.5, inplace=False)
#     )
#     (dense_1): Sequential(
#       (projection): Linear(in_features=200, out_features=80, bias=True)
#       (batch_norm): Identity()
#       (act_fn): ReLU()
#       (dropout): Dropout(p=0.5, inplace=False)
#     )
#   )
#   (final_dense): Sequential(
#     (0): Linear(in_features=80, out_features=1, bias=True)
#     (1): Sigmoid()
#   )
# )

model._associate_language(smiles_language)
model._associate_language(protein_language)

model_path = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/tutorial/TITAN_small/"
model_file = os.path.join(
    model_path, 'weights', 'best_ROC-AUC_bimodal_mca.pt'
)

if os.path.isfile(model_file):
    print('Found existing model, restoring now...')
    model.load(model_file, map_location=device)
    print(f'model loaded: {model_file}')


# Measure validation performance
loss_validation = []
# Sets the model to evaluation mode. This disables certain behaviors like dropout and batch normalization updates, ensuring consistent predictions during validation/testing.
model.eval() 
# torch.no_grad(): Disables gradient calculations, saving memory and computation since gradients are unnecessary during evaluation.
# test_loss = 0: A variable to accumulate the total loss across all batches in the test dataset.
# predictions: A list to store the model's output (predicted values) for each batch.
# labels: A list to store the actual labels (y) for comparison with predictions.
# tcr_attention: A list to collect attention scores (receptor-level attention) from the model for each batch, stored in pred_dict['receptor_attention'].
# tcr_id: A list to store receptor identifiers (or input receptor data) for each batch.

# model.loss: A loss function (e.g., Binary Cross Entropy Loss, BCELoss) defined in the model calculates the difference between the predicted (y_hat) and true labels (y).
# test_loss += loss.item(): The batch loss is added to test_loss, which accumulates the total validation loss.

with torch.no_grad():
    test_loss = 0
    predictions = []
    labels = []
    tcr_attention = []
    tcr_id = []
    for ind, (ligand, receptors, y) in enumerate(test_loader):
        y_hat, pred_dict = model(ligand.to(device), receptors.to(device))
        predictions.append(y_hat)
        labels.append(y.clone())
        loss = model.loss(y_hat, y.to(device))
        test_loss += loss.item()
        tcr_attention.append(pred_dict['receptor_attention'])
        tcr_id.append(receptors)

predictions = torch.cat(predictions, dim=0).flatten().cpu().numpy()
labels = torch.cat(labels, dim=0).flatten().cpu().numpy()
loss_validation.append(test_loss / len(test_loader))
tcr_attention = torch.cat(tcr_attention, dim=0).cpu().numpy()
tcr_id = torch.cat(tcr_id, dim=0).cpu().numpy()

test_loss = test_loss / len(test_loader)
fpr, tpr, _ = roc_curve(labels, predictions)
precision, recall, _ = precision_recall_curve(labels, predictions)
test_roc_auc = auc(fpr, tpr)
print('ROC AUC', test_roc_auc)

# Plot Evaluation Results
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
ax1.plot(fpr, tpr)
ax1.set_title('ROC Curve')
ax1.set_ylabel('TPR')
ax1.set_xlabel('FPR')
ax2.plot(precision, recall)
ax2.set_title('PR Curve')
ax2.set_ylabel('Recall')
ax2.set_xlabel('Precision')
plt.tight_layout()
plt.savefig("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/tutorial/evaluation_result.pdf")


def find_colors(attention, text):
    normalize = mpl.colors.Normalize(
            vmin=min(attention),
            vmax=max(attention)
        )
    color_mapper = cm.ScalarMappable(
        norm=normalize, cmap=plt.cm.Oranges
    )
    colors = color_mapper.to_rgba(attention, alpha=1.)
    if text:
        colors = [(x[0]*255,x[1]*255,x[2]*255) for x in colors]
    return(colors)

from IPython.core.display import display, HTML
import html 

def html_escape(text):
    return html.escape(text)

def highlight_text(text, colors, start_underline, end_underline):
    highlighted_text = []
    for i in range(len(text)):
        if i in start_underline:
            highlighted_text.append('<u>'
                                    +'<span style="font-family:courier; background-color:rgba(' 
                                    + str(int(colors[i][0])) + ','
                                    + str(int(colors[i][1])) + ','
                                    + str(int(colors[i][2])) + ');">' 
                                    + html_escape(text[i])
                                    + '</span>')
        elif i in end_underline:
            highlighted_text.append('<span style="font-family:courier; background-color:rgba(' 
                                    + str(int(colors[i][0])) + ','
                                    + str(int(colors[i][1])) + ','
                                    + str(int(colors[i][2])) + ');">' 
                                    + html_escape(text[i])
                                    + '</span>'
                                    + '</u>')
        else:
            highlighted_text.append('<span style="font-family:courier;background-color:rgba(' 
                                + str(int(colors[i][0])) + ','
                                + str(int(colors[i][1])) + ','
                                + str(int(colors[i][2])) + ');">' 
                                + html_escape(text[i])
                                + '</span>')
    highlighted_text = ' '.join(highlighted_text)
    return highlighted_text

show_aas = 50
cutoff = 130-show_aas

# Extract Attention Scores
tokens = []
for i in range(len(tcr_id)):
    t = []
    for i,x in enumerate(tcr_id[i]):
        if i > cutoff:
            if protein_language.index_to_token[x] == '<PAD>':
                token = '0'
            elif protein_language.index_to_token[x] == '<START>' or protein_language.index_to_token[x] == '<STOP>':
                token = '1'
            else:
                token = protein_language.index_to_token[x]
            t.append(token)
    tokens.append(t)


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

# Create a PDF to save the sequences
with PdfPages('highlighted_sequences.pdf') as pdf:
    # Create a figure for the plot
    fig, ax = plt.subplots(figsize=(10, 2))  # Adjust the figure size as needed
    for i in [0, 1, 2, 3, 4, 5, 6, 7]:
        seq = ''.join(tokens[i])  # Join the tokens to form the sequence
        att = tcr_attention[i][cutoff:]  # Get attention values
        col = find_colors(att, True)  # Get colors based on attention
        # Generate the highlighted sequence as a list of colored regions
        colored_sequence = []
        for j, char in enumerate(seq):
            colored_sequence.append((char, col[j]))  # Associate each character with its corresponding color
        # Create a plot to show the colored sequence
        ax.clear()  # Clear the previous plot
        ax.set_xlim(0, len(seq))  # Set the x-axis limits based on the sequence length
        ax.set_ylim(0, 1)  # Set the y-axis limits (fixed as we're just displaying text)
        # Plot each character with its corresponding color
        for j, (char, color) in enumerate(colored_sequence):
            ax.text(j, 0.5, char, ha='center', va='center', fontsize=12, color=color)
        ax.axis('off')  # Turn off axis to only display the sequence
        ax.set_title(f'Sequence {i + 1}')  # Title with sequence index
        # Save the figure to the PDF
        pdf.savefig(fig)
        plt.close(fig)  # Close the figure to prepare for the next one
