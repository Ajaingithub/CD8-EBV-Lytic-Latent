from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
import pandas as pd

# Define function to convert peptide sequence to SMILES
def peptide_to_smiles(sequence):
    mol = Chem.MolFromFASTA(sequence)
    if mol:
        return Chem.MolToSmiles(mol)
    return None

# List of epitopes
epitopes = ["ALWALPHAA", "FMVFLQTHI", "GLCTLVAML", "LLDFVRFMGV", "YVLDHLIVV"]

# Convert each epitope to SMILES
smiles_list = {epitope: peptide_to_smiles(epitope) for epitope in epitopes}

# Print results
for epitope, smiles in smiles_list.items():
    print(f"{epitope}: {smiles}")

smiles_df = pd.DataFrame(smiles_list.values())
smiles_df.to_csv("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/epitope.smi",index=False)


