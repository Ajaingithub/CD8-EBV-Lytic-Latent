#### We are developing our model to identify TCR and epitope prediction 
### Extracting out the data from seurat object
'''
CD8_Ag_metadata <- CD8_Ag@meta.data[,c("TRA1_TRA2_TRB1_TRB2_cdraa","TRA1_or_TRA2_TRB1_no_TRB2","TRA1_TRA2_TRB1_TRB2_fwr_aa","EBV_BALF4",
                                       "EBV_BMLF1","EBV_BMRF1","EBV_BRLF1","EBV_EBNA1","EBV_EBNA3C","EBV_LMP2","VZV_IE62","VZV_IE63","Ag_range_10","max_Ag")]
savedir = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/"
write.table(CD8_Ag_metadata, paste0(savedir,"data/TCR_EBV_peptide.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
'''

### Performing in Python
# conda activate af_py310
import pandas as pd
savedir = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/"
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
peptide_numeric["max_Ag"] = df_cleaned["max_Ag"]

TRA_TRB=df_cleaned.loc[:,['TRA1_CDR3','TRB1_CDR3']]
TCR_peptide = pd.concat([TRA_TRB, peptide_numeric], axis=1)

### Removing out the TCR which has low probability for any of the peptide
peptide_columns = [col for col in TCR_peptide.columns if '_prob' in col]
TCR_peptide_ge30 = TCR_peptide[TCR_peptide[peptide_columns].ge(30).any(axis=1)]
print(df_cleaned.shape,TCR_peptide_ge30.shape)

TCR_peptide_ge30.to_csv('TCR_EBV_peptide_processed.csv', index=False)

TCR_peptide_ge30_clean = TCR_peptide_ge30.drop_duplicates(subset=["TRA1_CDR3", "TRB1_CDR3"])
TCR_peptide_ge30_clean.to_csv('TCR_EBV_peptide_processed_nodups.csv', index=False)

df['max_Ag_encoded'].to_csv

##### Running a model on the unique TCR data
import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import LabelEncoder

os.chdir("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/")
df = pd.read_csv("TCR_EBV_peptide_processed_nodups.csv")

# Encode the antigen names (target variable)
label_encoder = LabelEncoder()
df['max_Ag_encoded'] = label_encoder.fit_transform(df['max_Ag'])


def one_hot_encode_multiple_sequences(sequences):
    """
    Encodes multiple amino acid sequences as one-hot encoded arrays, handling sequences of different lengths.
    Args:
    sequences (list of str): List of amino acid sequences. Each sequence is a string of amino acid characters.
    Returns:
    np.array: 3D NumPy array with shape (num_sequences, max_sequence_length, num_amino_acids)
    """
    # Determine the maximum sequence length
    max_sequence_length = max(len(seq) for seq in sequences)
    # Initialize a 3D array with zeros for the one-hot encoding
    num_sequences = len(sequences)
    num_amino_acids = len(amino_acids)
    one_hot_encoded = np.zeros((num_sequences, max_sequence_length, num_amino_acids), dtype=int)
    # Loop through each sequence and encode
    for seq_idx, sequence in enumerate(sequences):
        for aa_idx, aa in enumerate(sequence):
            # Set the corresponding index for the amino acid in the one-hot encoding
            one_hot_encoded[seq_idx, aa_idx, aa_to_index[aa]] = 1
    return one_hot_encoded

TRA1_CDR3_encoded = one_hot_encode_multiple_sequences(df["TRA1_CDR3"])
TRB1_CDR3_encoded = one_hot_encode_multiple_sequences(df["TRB1_CDR3"])
TRA1_TRB1_combined = np.concatenate((TRA1_CDR3_encoded, TRB1_CDR3_encoded), axis=1)

# Check the shape of the resulting array
print(TRA1_TRB1_combined.shape)  

# Step 6: Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(TRA1_TRB1_combined, df['max_Ag_encoded'], test_size=0.2, random_state=42)
print(X_train.shape, X_test.shape, y_train.shape,y_test.shape)

# Step 7: Train the Random Forest classifier
clf_model = RandomForestClassifier(n_estimators=100, random_state=42)
clf_model.fit(X_train, y_train)

# Step 8: Make predictions
y_pred = clf_model.predict(X_test)

### SInce we have to flat which would loose the information we are using LSTM
np.save('TRA1_TRB1_combined.npy', TRA1_TRB1_combined)
np.save('TRA1_CDR3_encoded.npy', TRA1_CDR3_encoded)
np.save('TRB1_CDR3_encoded.npy', TRB1_CDR3_encoded)
df['max_Ag_encoded'].to_csv('max_Ag_encoded.csv', index=False)

### RUnning it on GPU gpudev1
# conda activate alphafold2
# Define LSTM model
import os
import numpy as np
import pandas as pd
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

os.chdir("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data")
TRA1_TRB1_combined = np.load('TRA1_TRB1_combined.npy')
TRA1_CDR3_encoded = np.load('TRA1_CDR3_encoded.npy')
TRB1_CDR3_encoded = np.load('TRB1_CDR3_encoded.npy')
max_Ag_encoded = pd.read_csv("max_Ag_encoded.csv", sep="\t")


#  Example 3D array (n_samples, n_timestamps, n_features)
# n_samples = 10
# n_timestamps = 5
# n_features = 3
# X = np.random.rand(n_samples, n_timestamps, n_features)  # Shape: (10, 5, 3)

# # Target variable (for simplicity, using random integers)
# y = np.random.choice([0, 1], size=n_samples)

# Split data into training and testing sets (80% training, 20% testing)
X_train, X_test, y_train, y_test = train_test_split(TRB1_CDR3_encoded, max_Ag_encoded, test_size=0.2, random_state=42)
print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)

model = Sequential()
model.add(LSTM(units=50, activation='relu', input_shape=(25, 20)))
model.add(Dense(1, activation='sigmoid'))  # Binary classification (0 or 1)

# Compile and train the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.fit(X_train, y_train, epochs=20, batch_size=2, verbose=1)

# Evaluate the model
y_pred = model.predict(X_test)
y_pred = (y_pred > 0.5).astype(int)

# Evaluate accuracy
print(f"Accuracy: {accuracy_score(y_test, y_pred)}")
# Accuracy: 0.002237970906378217


####### CNN
# Define CNN model
import numpy as np
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

model = Sequential()
model.add(Conv1D(filters=64, kernel_size=3, activation='relu', input_shape=(25, 20)))

# Optionally, add a MaxPooling1D layer to reduce dimensionality
model.add(MaxPooling1D(pool_size=2))
model.add(Flatten())
model.add(Dense(32, activation='relu'))
model.add(Dense(1, activation='sigmoid'))  # Binary classification (0 or 1)

# Compile and train the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.fit(X_train, y_train, epochs=20, batch_size=2, verbose=1)

# Evaluate the model
y_pred = model.predict(X_test)
y_pred = (y_pred > 0.5).astype(int)  # Convert probabilities to binary labels

# Evaluate accuracy
print(f"Accuracy: {accuracy_score(y_test, y_pred)}")