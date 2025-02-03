import random
import pandas as pd

# List of tuples (word, number)
data = pd.read_csv("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/TCR_epitope.txt", header=None)

data[['TCR', 'epitope']] = data[0].str.split('_', expand=True)
data.drop(columns=[0], inplace=True)
print(data)

words = list(data['TCR'])  # List of words
numbers = list(data['epitope'])  # List of numbers

# Shuffle the numbers
shuffled_numbers = numbers[:]
random.shuffle(shuffled_numbers)

# Ensure that no word gets its original number
while any(numbers[i] == shuffled_numbers[i] for i in range(len(numbers))):
    random.shuffle(shuffled_numbers)

# Pair words with their new shuffled numbers
shuffled_data = list(zip(words, shuffled_numbers))

# Print the shuffled result
for word, num in shuffled_data:
    print(f"{word} {num}")

data = pd.read_csv("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/full_titan_dataset.txt", header = None, sep = "\t")
data.columns = ["epitopes","TCRs","labels"]

# Shuffle the 'epitopes' and 'TCRs' columns
shuffled_epitopes = data['epitopes'].sample(frac=1).reset_index(drop=True)
shuffled_TCRs = data['TCRs'].sample(frac=1).reset_index(drop=True)


# Create a new dataframe with shuffled values and labels set to 0
shuffled_data = pd.DataFrame({
    'epitopes': shuffled_epitopes,
    'TCRs': shuffled_TCRs,
    'labels': [0] * len(data)
})

merged_data = pd.merge(data, shuffled_data, on=['epitopes', 'TCRs'], how='inner')

merged_data.to_csv("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/same_TCR_epitope.csv", index=False)
# Remove duplicates from original and shuffled data
original_no_duplicates = data[~data.isin(merged_data)].dropna()
shuffled_no_duplicates = shuffled_data[~shuffled_data.isin(merged_data)].dropna()
