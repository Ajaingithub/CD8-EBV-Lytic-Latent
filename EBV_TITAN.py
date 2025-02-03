#### Running a TITAN for our constructed EBV dataset

python3 /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/scripts/flexible_training.py \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/Ag_train.csv \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/Ag_test.csv  \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/tcr.csv \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/epitope.smi \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/trained_model \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/tutorial/data/params_small.json \
All_Ag \
bimodal_mca

# INFO:VZV_IE62:	 **** TRAINING ****   Epoch [50/100], loss: 0.00086. This took 6.5 secs.
# INFO:VZV_IE62:	 **** TESTING **** Epoch [50/100], loss: 0.00011, ROC-AUC: 1.000, Average precision: 1.000.
# INFO:VZV_IE62:	 New best performance in "loss" with value : 0.0001058 in epoch: 49
# INFO:VZV_IE62:== Epoch [50/100] ==

# python3 /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/scripts/flexible_training.py \
# /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/trained_model_train.csv \
# /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/trained_model_test.csv \
# /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/tcr.csv \
# /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/epitopes.smi \
# trained_model \
# /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/tutorial/data/params_small.json \
# full_dataset \
# bimodal_mca


python3 /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/scripts/flexible_model_eval.py \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/Ag_train.csv \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/tcr.csv \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/epitope.smi \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/trained_model/full_dataset \
bimodal_mca \
Ag_evaluate

python3 /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/scripts/flexible_model_eval.py \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/trained_model_test.csv \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/tcr.csv \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/epitopes.smi \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/trained_model/full_dataset \
bimodal_mca \
tool_test_eval

python3 /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/scripts/flexible_model_eval.py \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/shuffled_test.txt \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/tcr.csv \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/datasets/epitopes.smi \
/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/trained_model/full_dataset \
bimodal_mca \
tool_test_eval

#### EBV TITAN training
### This will remove all the variables except some buildins
for var in list(globals().keys()):
    if var not in ["__builtins__", "__file__", "__name__", "__doc__", "__package__"]:
        del globals()[var]

"""Train Affinity predictor model."""
import argparse
import json
import logging
import os
import sys
from time import time

import numpy as np
import torch
from paccmann_predictor.models import MODEL_FACTORY
from paccmann_predictor.utils.hyperparams import OPTIMIZER_FACTORY
from paccmann_predictor.utils.utils import get_device
from pytoda.datasets import (
    DrugAffinityDataset, ProteinProteinInteractionDataset
)
from pytoda.proteins import ProteinFeatureLanguage, ProteinLanguage
from pytoda.smiles.smiles_language import SMILESTokenizer
from sklearn.metrics import (
    auc, average_precision_score, precision_recall_curve, roc_curve
)
from pytoda.smiles import metadata

torch.manual_seed(123456)

# setup logging
logging.basicConfig(stream=sys.stdout)

data_dir = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/"
train_affinity_filepath = os.path.join(data_dir, "train.csv")
test_affinity_filepath = os.path.join(data_dir, "test.csv")
receptor_filepath = os.path.join(data_dir, "tcr.csv")
ligand_filepath = os.path.join(data_dir, "epitope.smi")
model_path = os.path.join(data_dir, "EBV_trained_model")
params_filepath = os.path.join("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/TITAN/tutorial/data/params_small.json")
training_name = "EBV_CDR3"
model_type = "bimodal_mca"

# yapf: enable

# def main(
#     train_affinity_filepath, test_affinity_filepath, receptor_filepath,
#     ligand_filepath, model_path, params_filepath, training_name, model_type
# ):

logger = logging.getLogger(f'{training_name}')
logger.setLevel(logging.DEBUG)
# Process parameter file:
params = {}
with open(params_filepath) as fp:
    params.update(json.load(fp))

    # Create model directory and dump files
model_dir = os.path.join(model_path, training_name)
os.makedirs(os.path.join(model_dir, 'weights'), exist_ok=True)
os.makedirs(os.path.join(model_dir, 'results'), exist_ok=True)
with open(os.path.join(model_dir, 'model_params.json'), 'w') as fp:
    json.dump(params, fp, indent=4)

device = get_device()
print(device)
# Load languages
smiles_language_filepath = os.path.join(
    os.sep,
    *metadata.__file__.split(os.sep)[:-1], 'tokenizer'
    )
smiles_language = SMILESTokenizer.from_pretrained(smiles_language_filepath)
smiles_language.set_encoding_transforms(
    randomize=None,
    add_start_and_stop=params.get('ligand_start_stop_token', True),
    padding=params.get('ligand_padding', True),
    padding_length=params.get('ligand_padding_length', True),
    )
smiles_language.set_smiles_transforms(
    augment=params.get('augment_smiles', False),
    canonical=params.get('smiles_canonical', False),
    kekulize=params.get('smiles_kekulize', False),
    all_bonds_explicit=params.get('smiles_bonds_explicit', False),
    all_hs_explicit=params.get('smiles_all_hs_explicit', False),
    remove_bonddir=params.get('smiles_remove_bonddir', False),
    remove_chirality=params.get('smiles_remove_chirality', False),
    selfies=params.get('selfies', False),
    sanitize=params.get('sanitize', False)
    )

if params.get('receptor_embedding', 'learned') == 'predefined':
    protein_language = ProteinFeatureLanguage(
        features=params.get('predefined_embedding', 'blosum')
        )
else:
    protein_language = ProteinLanguage()

if params.get('ligand_embedding', 'learned') == 'one_hot':
    logger.warning(
        'ligand_embedding_size parameter in param file is ignored in '
        'one_hot embedding setting, ligand_vocabulary_size used instead.'
        )
if params.get('receptor_embedding', 'learned') == 'one_hot':
    logger.warning(
        'receptor_embedding_size parameter in param file is ignored in '
        'one_hot embedding setting, receptor_vocabulary_size used instead.'
        )

# Prepare the dataset
logger.info("Start data preprocessing...")

# Check if peptide as SMILES or as aa
pepname, pep_extension = os.path.splitext(ligand_filepath)
if pep_extension == '.csv':
    logger.info(
        'Ligand file has extension .csv \n'
        'Please make sure ligand is provided as amino acid sequence.'
    )
    train_dataset = ProteinProteinInteractionDataset(
        sequence_filepaths=[[ligand_filepath], [receptor_filepath]],
        entity_names=['ligand_name', 'sequence_id'],
    labels_filepath=train_affinity_filepath,
    annotations_column_names=['label'],
    protein_languages=protein_language,
    padding_lengths=[
        params.get('ligand_padding_length', None),
        params.get('receptor_padding_length', None)
    ],
    paddings=params.get('ligand_padding', True),
    add_start_and_stops=params.get('add_start_stop_token', True),
    augment_by_reverts=params.get('augment_protein', False),
    randomizes=params.get('randomize', False),
    iterate_datasets=True
    )
    train_loader = torch.utils.data.DataLoader(
        dataset=train_dataset,
        batch_size=params['batch_size'],
        shuffle=True,
        drop_last=True,
        num_workers=params.get('num_workers', 0)
        )

    test_dataset = ProteinProteinInteractionDataset(
        sequence_filepaths=[[ligand_filepath], [receptor_filepath]],
        entity_names=['ligand_name', 'sequence_id'],
        labels_filepath=test_affinity_filepath,
        annotations_column_names=['label'],
        protein_languages=protein_language,
        padding_lengths=[
            params.get('ligand_padding_length', None),
            params.get('receptor_padding_length', None)
        ],
        paddings=params.get('ligand_padding', True),
        add_start_and_stops=params.get('add_start_stop_token', True),
        augment_by_reverts=params.get('augment_test_data', False),
        randomizes=False,
        iterate_datasets=True
    )

    test_loader = torch.utils.data.DataLoader(
        dataset=test_dataset,
        batch_size=params['batch_size'],
        shuffle=True,
        drop_last=True,
        num_workers=params.get('num_workers', 0)
    )

    params.update({
        'ligand_vocabulary_size': protein_language.number_of_tokens,
        'receptor_vocabulary_size': protein_language.number_of_tokens,
        'ligand_as': 'amino acids'
    })  # yapf: disable

    logger.info(
        f'ligand_vocabulary_size {protein_language.number_of_tokens}, '
        f'receptor_vocabulary_size {protein_language.number_of_tokens}'
    )

    logger.info(
        f'Training dataset has {len(train_dataset)} samples, test set has '
        f'{len(test_dataset)}.'
    )

elif pep_extension == '.smi':
    logger.info(
        'Ligand file has extension .smi \n'
        'Please make sure ligand is provided as SMILES.'
    )

    print(protein_language.max_token_sequence_length)  # Keep this print for debugging if needed
        
    # Assemble datasets
    train_dataset = DrugAffinityDataset(
        drug_affinity_filepath=train_affinity_filepath,
        smi_filepath=ligand_filepath,
        protein_filepath=receptor_filepath,
        smiles_language=smiles_language,
        protein_language=protein_language,
        smiles_padding=params.get('ligand_padding', True),
        smiles_padding_length=params.get('ligand_padding_length', None),
        smiles_add_start_and_stop=params.get('ligand_add_start_stop', True),
        smiles_augment=params.get('augment_smiles', False),
        smiles_canonical=params.get('smiles_canonical', False),
        smiles_kekulize=params.get('smiles_kekulize', False),
        smiles_all_bonds_explicit=params.get('smiles_bonds_explicit', False),
        smiles_all_hs_explicit=params.get('smiles_all_hs_explicit', False),
        smiles_remove_bonddir=params.get('smiles_remove_bonddir', False),
        smiles_remove_chirality=params.get('smiles_remove_chirality', False),
        smiles_selfies=params.get('selfies', False),
        protein_amino_acid_dict=params.get('protein_amino_acid_dict', 'iupac'),
        protein_padding=params.get('receptor_padding', True),
        protein_padding_length=params.get('receptor_padding_length', None),
        protein_add_start_and_stop=params.get('receptor_add_start_stop', True),
        protein_augment_by_revert=params.get('augment_protein', False),
        drug_affinity_dtype=torch.float,
        backend='eager',
        iterate_dataset=True
    )

    train_loader = torch.utils.data.DataLoader(
        dataset=train_dataset,
        batch_size=params['batch_size'],
        shuffle=True,
        drop_last=True,
        num_workers=params.get('num_workers', 0)
    )

    # test_dataset = DrugAffinityDataset(
    #     drug_affinity_filepath=test_affinity_filepath,
    #     smi_filepath=ligand_filepath,
    #     protein_filepath=receptor_filepath,
    #     smiles_language=smiles_language,
    #     protein_language=protein_language,
    #     smiles_padding=params.get('ligand_padding', True),
    #     smiles_padding_length=params.get('ligand_padding_length', None),
    #     smiles_add_start_and_stop=params.get('ligand_add_start_stop', True),
    #     smiles_augment=False,  # Augmentation should usually be off for test data
    #     smiles_canonical=params.get('test_smiles_canonical', False),
    #     smiles_kekulize=params.get('smiles_kekulize', False),
    #     smiles_all_bonds_explicit=params.get('smiles_bonds_explicit', False),
    #     smiles_all_hs_explicit=params.get('smiles_all_hs_explicit', False),
    #     smiles_remove_bonddir=params.get('smiles_remove_bonddir', False),
    #     smiles_remove_chirality=params.get('smiles_remove_chirality', False),
    #     smiles_selfies=params.get('selfies', False),
    #     protein_amino_acid_dict=params.get('protein_amino_acid_dict', 'iupac'),
    #     protein_padding=params.get('receptor_padding', True),
    #     protein_padding_length=params.get('receptor_padding_length', None),
    #     protein_add_start_and_stop=params.get('receptor_add_start_stop', True),
    #     protein_augment_by_revert=False,  # Augmentation should usually be off for test data
    #     drug_affinity_dtype=torch.float,
    #     backend='eager',
    #     iterate_dataset=True
    # )

    test_dataset = DrugAffinityDataset(
        drug_affinity_filepath=test_affinity_filepath,
        smi_filepath=ligand_filepath,
        protein_filepath=receptor_filepath,
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

    logger.info(
        f'Training dataset has {len(train_dataset)} samples, test set has '
        f'{len(test_dataset)}.'
    )

    test_loader = torch.utils.data.DataLoader(
        dataset=test_dataset,
        batch_size=params['batch_size'],
        shuffle=False,  # Shuffle test data if needed
        drop_last=False,
        num_workers=params.get('num_workers', 0)
    )

    params.update({
        'ligand_vocabulary_size': smiles_language.number_of_tokens,
        'receptor_vocabulary_size': protein_language.number_of_tokens,
        'ligand_as': 'smiles'
    })  # yapf: disable

    logger.info(
        f'ligand_vocabulary_size {smiles_language.number_of_tokens}, '
        f'receptor_vocabulary_size {protein_language.number_of_tokens}'
    )

else:
    raise ValueError(
        f"Choose pep_filepath with extension .csv or .smi, given was {pep_extension}"
    )
    
save_top_model = os.path.join(model_dir, 'weights/{}_{}_{}.pt')

model_fn = params.get('model_fn', model_type)
model = MODEL_FACTORY[model_fn](params).to(device)
model._associate_language(smiles_language)
model._associate_language(protein_language)

smiles_language.save_pretrained(model_dir)
protein_language.save(os.path.join(model_dir, 'protein_language.pkl'))

# Define optimizer
min_loss, max_roc_auc = 100, 0
optimizer = OPTIMIZER_FACTORY[params.get('optimizer', 'adam')](
    model.parameters(),
    lr=params.get('lr', 0.001),
    weight_decay=params.get('weight_decay', 0.001)
)

num_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
params.update({'number_of_parameters': num_params})
logger.info(f'Number of parameters: {num_params}')
logger.info(f'Model: {model}')

# Overwrite params.json file with updated parameters.
with open(os.path.join(model_dir, 'model_params.json'), 'w') as fp:
    json.dump(params, fp)

# Start training
logger.info('Training about to start...\n')
t = time()
loss_training = []
loss_validation = []

model.save(save_top_model.format('epoch', '0', model_fn))

for epoch in range(params['epochs']):
    model.train()
    logger.info(f"== Epoch [{epoch+1}/{params['epochs']}] ==")  # Corrected epoch numbering in log
    train_loss = 0

    for ind, (ligand, receptors, y) in enumerate(train_loader):
        torch.cuda.empty_cache()
        if ind % 10 == 0:
            logger.info(f'Batch {ind}/{len(train_loader)}')
        y_hat, pred_dict = model(ligand.to(device), receptors.to(device))
        loss = model.loss(y_hat, y.to(device))
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        train_loss += loss.item()

    logger.info(
        "\t **** TRAINING ****   "
        f"Epoch [{epoch+1}/{params['epochs']}], "  # Corrected epoch numbering in log
        f"loss: {train_loss / len(train_loader):.5f}. "
        f"This took {time() - t:.1f} secs."
        )

    t = time()

    # Measure validation performance
    model.eval()
    with torch.no_grad():
        test_loss = 0
        predictions = []
        labels = []
        for ind, (ligand, receptors, y) in enumerate(test_loader):
            torch.cuda.empty_cache()
            y_hat, pred_dict = model(ligand.to(device), receptors.to(device))
            predictions.append(y_hat)
            labels.append(y.clone())
            loss = model.loss(y_hat, y.to(device))
            test_loss += loss.item()
            # predictions.append(y_hat.cpu().detach())
            # predictions.append(y_hat.cpu().detach().numpy())  # Converts to NumPy for later analysis
    
    print(predictions)
    # logger.info("Predictions:", predictions)  # Keep this for debugging if needed.  Consider logging instead of print.
    
    predictions = torch.cat(predictions, dim=0).flatten().cpu().numpy()
    labels = torch.cat(labels, dim=0).flatten().cpu().numpy()
    loss_validation.append(test_loss / len(test_loader))
    loss_training.append(train_loss / len(train_loader))

    test_loss = test_loss / len(test_loader)
    fpr, tpr, _ = roc_curve(labels, predictions)
    test_roc_auc = auc(fpr, tpr)

    # calculations for visualization plot
    precision, recall, _ = precision_recall_curve(labels, predictions)
    avg_precision = average_precision_score(labels, predictions)
    
    logger.info(
        f"\t **** TESTING **** Epoch [{epoch+1}/{params['epochs']}], "  # Corrected epoch numbering in log
        f"loss: {test_loss:.5f}, ROC-AUC: {test_roc_auc:.3f}, "
        f"Average precision: {avg_precision:.3f}."
    )

    def save(path, metric, typ, val=None):
        model.save(path.format(typ, metric, model_fn))
        info = {
            'best_roc_auc': str(max_roc_auc),
            'test_loss': str(min_loss)
        }
        with open(
            os.path.join(model_dir, 'results', metric + '.json'), 'w'
        ) as f:
            json.dump(info, f)
        np.save(
            os.path.join(model_dir, 'results', metric + '_preds.npy'),
            np.vstack([predictions, labels])
        )
        if typ == 'best':
            logger.info(
                f'\t New best performance in "{metric}"'
                f' with value : {val:.7f} in epoch: {epoch}'
            )

    if test_roc_auc > max_roc_auc:
        max_roc_auc = test_roc_auc
        save(save_top_model, 'ROC-AUC', 'best', max_roc_auc)
        ep_roc = epoch
        roc_auc_loss = test_loss
        roc_auc_pr = avg_precision

    if test_loss < min_loss:
        min_loss = test_loss
        save(save_top_model, 'loss', 'best', min_loss)
        ep_loss = epoch
        loss_roc_auc = test_roc_auc
    if (epoch + 1) % params.get('save_model', 100) == 0:
        save(save_top_model, 'epoch', str(epoch))

logger.info(
    'Overall best performances are: \n \t'
    f'Loss = {min_loss:.4f} in epoch {ep_loss} '
    f'\t (ROC-AUC was {loss_roc_auc:4f}) \n \t'
    f'ROC-AUC = {max_roc_auc:.4f} in epoch {ep_roc} '
    f'\t (Loss was {roc_auc_loss:4f})'
)

save(save_top_model, 'training', 'done')
logger.info('Done with training, models saved, shutting down.')

np.save(
    os.path.join(model_dir, 'results', 'loss_training.npy'), loss_training
)
np.save(
    os.path.join(model_dir, 'results', 'loss_validation.npy'),
    loss_validation
)

    # save best results
result_file = os.path.join(model_path, 'results_overview.csv')
with open(result_file, "a") as myfile:
    myfile.write(
        f'{training_name},{max_roc_auc:.4f},{roc_auc_pr:.4f},{roc_auc_loss:.4f},{ep_roc} \n \t'
    )