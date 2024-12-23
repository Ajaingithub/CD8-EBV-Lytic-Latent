#### We are developing our model to identify TCR and epitope prediction 
### Extracting out the data from seurat object
'''
CD8_Ag_metadata <- CD8_Ag@meta.data[,c("TRA1_TRA2_TRB1_TRB2_cdraa","TRA1_or_TRA2_TRB1_no_TRB2","TRA1_TRA2_TRB1_TRB2_fwr_aa","EBV_BALF4",
                                       "EBV_BMLF1","EBV_BMRF1","EBV_BRLF1","EBV_EBNA1","EBV_EBNA3C","EBV_LMP2","VZV_IE62","VZV_IE63","Ag_range_10","max_Ag")]
savedir = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/"
write.table(CD8_Ag_metadata, paste0(savedir,"data/TCR_EBV_peptide.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
'''

### Performing in Python
import pandas as pd
df = pd.read_csv(f"{savedir}/data/TCR_EBV_peptide.txt", sep='\t')
df_cleaned = df.dropna() ### Removed Nan
print(df.shape, df_cleaned.shape)

### Separating out TRA and TRB based on --
df_cleaned[["TRA_aa","TRB_aa"]] = df_cleaned['TRA1_TRA2_TRB1_TRB2_cdraa'].str.split('--', expand=True)

### Keeping only TRA1 and TRB1
df_cleaned['TRA1_aa'] = df_cleaned['TRA_aa'].str.replace(r'-.*', '', regex=True) 
df_cleaned['TRB1_aa'] = df_cleaned['TRB_aa'].str.replace(r'-.*', '', regex=True)

### Keeping only TRA1 and TRB1 CDR3
# ^.*?_.*?_: This matches everything from the start of the string (^), followed by any characters up to the first underscore (.*?_),
# then everything up to the second underscore (.*?_), and removes that part of the string.
df_cleaned['TRA1_CDR3'] = df_cleaned['TRA1_aa'].str.replace(r'^.*?_.*?_', '', regex=True) 
df_cleaned['TRB1_CDR3'] = df_cleaned['TRB1_aa'].str.replace(r'^.*?_.*?_', '', regex=True)

# Identify columns starting with 'EBV' or 'VZV'
ebv_vzv_columns = [col for col in df_cleaned.columns if col.startswith('EBV') or col.startswith('VZV')]

# Loop through each 'EBV'/'VZV' column to extract the last value
for col in ebv_vzv_columns:
    # Extract the last part after the last hyphen from each value in the column
    new_col_name = f'{col}_prob'  # Creating the new column name (e.g., EBV_BALF4_prob)
    df_cleaned[new_col_name] = df_cleaned[col].str.split('-').str[-1]

required_col = [col for col in df_cleaned.columns if col.endswith('_prob')]
peptide = df_cleaned.loc[:,required_col]
peptide_numeric = peptide.apply(pd.to_numeric, errors='coerce')
peptide_numeric["Ag_range_10"] = df_cleaned["Ag_range_10"]
peptide_numeric["Ag_range_10"] = df_cleaned["Ag_range_10"]


TRA_TRB=df_cleaned.loc[:,['TRA1_CDR3','TRB1_CDR3']]
TCR_peptide = pd.concat([TRA_TRB, peptide_numeric], axis=1)

### Removing out the TCR which has low probability for any of the peptide
peptide_columns = [col for col in TCR_peptide.columns if '_prob' in col]
TCR_peptide_ge30 = TCR_peptide[TCR_peptide[peptide_columns].ge(30).any(axis=1)]
print(df_cleaned.shape,TCR_peptide_ge30.shape)

TCR_peptide_ge30.to_csv('TCR_EBV_peptide_processed.csv', index=False)

