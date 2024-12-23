import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import LabelEncoder

# Step 1: Load the Data
sequences_file = 'TRA1_CDR3.txt'  # Path to the sequences file
sequences_df = pd.read_csv(sequences_file, header=None, names=['ID', 'Sequence'], delimiter='\t')
sequences_df = sequences_df.dropna(subset=['Sequence'])
#label_encoder = LabelEncoder()
#sequences_df['Encoded_Sequence'] = label_encoder.fit_transform(sequences_df['Sequence'])

def one_hot_encode(sequence):
    # Define the alphabet of amino acids
    aa_alphabet = list('ARNDCQEGHILKMFPSTWYV')
    encoding = {aa: i for i, aa in enumerate(aa_alphabet)}
    print(encoding)
    # Create a binary matrix for the sequence
    encoded_seq = np.zeros(len(aa_alphabet), dtype=int)
    for aa in sequence:
        if aa in encoding:
            encoded_seq[encoding[aa]] = 1
    return encoded_seq

# Apply one-hot encoding to each sequence
encoded_sequences = np.array([one_hot_encode(seq) for seq in sequences_df['Sequence']])

binding_scores_file = 'TCR_Peptide_value_acast.txt'  # Path to the binding scores file
binding_scores_df = pd.read_csv(binding_scores_file, sep='\t')
binding_scores_df = binding_scores_df[binding_scores_df.drop(columns=['ColNum']).gt(30).any(axis=1)]

# Step 3: Merge Data with Binding Scores
# Merge sequences and binding scores data based on the 'ID' column
merged_df = pd.merge(sequences_df, binding_scores_df, left_on='ID', right_on='ColNum')

binding_score_columns = binding_scores_df.columns[1:]
target_classes = merged_df[binding_score_columns].idxmax(axis=1)  # Get the column with the max value for each row
target_classes = target_classes.map(lambda x: binding_score_columns.get_loc(x))

X = merged_df['Encoded_Sequence'].values.reshape(-1, 1)  # Features: encoded sequence
y = target_classes  # Target: the index of the highest binding score

# Step 6: Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Step 7: Train the Random Forest classifier
clf_model = RandomForestClassifier(n_estimators=100, random_state=42)
clf_model.fit(X_train, y_train)

# Step 8: Make predictions
y_pred = clf_model.predict(X_test)

# Step 9: Evaluate the model
from sklearn.metrics import accuracy_score
accuracy = accuracy_score(y_test, y_pred)
print(f'Accuracy: {accuracy}')


# Step 7: Evaluate the Model
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"Mean Squared Error: {mse}")
print(f"R-squared: {r2}")

# Optional: Feature Importance
import matplotlib.pyplot as plt

# Plot feature importance
feature_importance = model.feature_importances_
plt.figure(figsize=(10, 6))
plt.bar(range(len(feature_importance)), feature_importance)
plt.title("Feature Importance")
plt.xlabel("Amino Acid Position")
plt.ylabel("Importance")
plt.savefig('feature_importance_plot.png')  
# plt.savefig('feature_importance_plot.pdf')  # Save as PDF
# plt.savefig('feature_importance_plot.svg')  # Save as SVG
plt.close()


def calculate_accuracy(y_true, y_pred, tolerance=0.1):
    # Calculate the absolute error
    absolute_error = np.abs(y_true - y_pred)
    
    # Calculate the percentage of predictions within the tolerance (e.g., 10%)
    accuracy_mask = absolute_error <= tolerance * y_true  # Tolerance threshold based on true value
    accuracy_percentage = np.mean(accuracy_mask) * 100  # Percentage of "accurate" predictions
    
    return accuracy_percentage

# Calculate the custom accuracy on the test set
accuracy_percentage = calculate_accuracy(y_test, y_pred, tolerance=0.1)

# Print the accuracy
print(f"Accuracy Percentage: {accuracy_percentage:.2f}%")

