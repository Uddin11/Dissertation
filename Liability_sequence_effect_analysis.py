import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# Define liability patterns
liability_patterns = {
    'Deamidation': re.compile(r'N(?=G)'),
    'Oxidation': re.compile(r'M|W'),
    'Glycosylation': re.compile(r'N[^P][ST]'),
    'Aggregation': re.compile(r'F{2,}|W{2,}|Y{2,}|L{2,}|V{2,}|I{2,}|M{2,}|A{2,}')
}

# Function to identify liability sequences in an antibody sequence
def identify_liabilities(sequence):
    liability_counts = {key: len(pattern.findall(sequence)) for key, pattern in liability_patterns.items()}
    return liability_counts

# Directories containing the ANARCI files
anarci_directories = [
    '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/bound_renumbered',
    '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/unbound_renumbered'
]

# List of specific antibodies to investigate
antibodies_of_interest = [
    '1a14', '1a3l', '1a3r', '1a4k', '1a6v', '1a6w', '1acy', '1adq', '1afv', '1ai1', 
    '1aj7', '1axs', '1baf', '1bfv', '1bgx', '1bj1', '1bql', '1c1e', '1c5c', '1cbv', 
    '1ce1', '1cf8', '1cfv', '1clz', '1cly', '1cu4', '1cz8', '1d6v', '1dbb', '1dbj', 
    '1dbk', '1dbm'
]

# Extract liabilities from sequences
liability_results = []
for anarci_dir in anarci_directories:
    for anarci_file in os.listdir(anarci_dir):
        if anarci_file.endswith('.renumbered'):
            antibody_name = anarci_file.split('.')[0]
            if antibody_name in antibodies_of_interest:
                file_path = os.path.join(anarci_dir, anarci_file)
                with open(file_path, 'r') as file:
                    sequence = ''
                    for line in file:
                        parts = line.strip().split()
                        if len(parts) > 2:  # must be atleast three in the line
                            sequence += parts[2]  # Extract the sequence part
                    liabilities = identify_liabilities(sequence)
                    total_liabilities = sum(liabilities.values())
                    liability_results.append({'Antibody': antibody_name, 'Total Liabilities': total_liabilities, **liabilities})

# Convert liability results to a DataFrame
liability_df = pd.DataFrame(liability_results)

# Read RMSD values from CSV files
rmsd_directory = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/'
rmsd_data = []
for rmsd_file in os.listdir(rmsd_directory):
    if rmsd_file.endswith('_rmsd_results.csv'):
        antibody_name = rmsd_file.split('_')[0]
        if antibody_name in antibodies_of_interest:
            df = pd.read_csv(os.path.join(rmsd_directory, rmsd_file))
            average_rmsd = df['RMSD Value'].mean()
            rmsd_data.append({'Antibody': antibody_name, 'Average RMSD': average_rmsd})

# Convert RMSD data to a DataFrame
rmsd_df = pd.DataFrame(rmsd_data)

# Merge liability data with RMSD data
if 'Antibody' in liability_df.columns and 'Antibody' in rmsd_df.columns:
    merged_df = pd.merge(liability_df, rmsd_df, on='Antibody')

    # Plot the amount of liability sequences
    plt.figure(figsize=(12, 6))
    plt.bar(merged_df['Antibody'], merged_df['Total Liabilities'], color='red', label='Total Liabilities')
    plt.xlabel('Antibody')
    plt.ylabel('Total Liability Sequences')
    plt.title('Total Liability Sequences in Antibody Sequences')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(os.path.join(rmsd_directory, 'liability_sequences_overview.png'))
    plt.show()

    # Plot the comparison between liability sequences and RMSD
    plt.figure(figsize=(12, 6))
    plt.scatter(merged_df['Total Liabilities'], merged_df['Average RMSD'], color='blue')
    plt.xlabel('Total Liability Sequences')
    plt.ylabel('Average RMSD')
    plt.title('Comparison of Total Liability Sequences and RMSD')
    plt.tight_layout()
    plt.savefig(os.path.join(rmsd_directory, 'liability_vs_rmsd.png'))
    plt.show()

    # Save merged data to CSV for further analysis
    merged_df.to_csv(os.path.join(rmsd_directory, 'merged_liability_rmsd_data.csv'), index=False)
else:
    print("Error: 'Antibody' column missing in one of the DataFrames.")




# Compare TAP analysis metric with RMSD
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Paths to the TAP analysis and RMSD data
tap_analysis_path = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/Liability_sequence/TAP_Analysis.csv'
rmsd_directory = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/'
output_directory = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/Liability_sequence'

# Load TAP analysis data
try:
    tap_df = pd.read_csv(tap_analysis_path)
    print("TAP data loaded successfully.")
except FileNotFoundError:
    print(f"Error: File not found at {tap_analysis_path}")
    raise

# Check if 'Value' column exists
if 'Value' not in tap_df.columns:
    raise ValueError("'Value' column is missing in the TAP data.")

# Pivot the data to wide format
try:
    tap_wide_df = tap_df.pivot(index='Name', columns='Metric', values='Value').reset_index()
    print("Pivot operation successful.")
except KeyError as e:
    print(f"Error during pivot operation: {e}")
    raise

# Load RMSD data
rmsd_data = []
for rmsd_file in os.listdir(rmsd_directory):
    if rmsd_file.endswith('_rmsd_results.csv'):
        antibody_name = rmsd_file.split('_')[0]
        df = pd.read_csv(os.path.join(rmsd_directory, rmsd_file))
        average_rmsd = df['RMSD Value'].mean()
        rmsd_data.append({'Name': antibody_name, 'Average RMSD': average_rmsd})

rmsd_df = pd.DataFrame(rmsd_data)

# Merge TAP metrics with RMSD data
try:
    merged_df = pd.merge(tap_wide_df, rmsd_df, on='Name')
    print("Data merged successfully.")
except KeyError as e:
    print(f"Error during merge operation: {e}")
    raise

# Ensure all relevant columns are numeric
metrics_columns = ['Total CDR Length', 'CDR Vicinity PSH Score (Kyte & Doolittle)', 'CDR Vicinity PPC Score', 'CDR Vicinity PNC Score', 'SFvCSP Score']
valid_metrics = []

for column in metrics_columns:
    if column in merged_df.columns:
        merged_df[column] = pd.to_numeric(merged_df[column], errors='coerce')
        valid_metrics.append(column)
    else:
        print(f"Column '{column}' not found in the merged DataFrame. Skipping this metric.")

# Drop rows where any valid metric column is still NaN
merged_df.dropna(subset=valid_metrics, inplace=True)

# Plot pairwise comparisons with line of best fit
for metric in valid_metrics:
    sns.lmplot(x=metric, y='Average RMSD', data=merged_df, aspect=2)
    plt.title(f'{metric} vs. RMSD with Line of Best Fit')
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, f'{metric}_vs_RMSD_with_fit.png'))
    plt.show()

# Analyze correlations between TAP metrics and RMSD
correlation_matrix = merged_df.corr()

# Plot correlation heatmap and save it
plt.figure(figsize=(12, 10))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f")
plt.title('Correlation between TAP Metrics and RMSD')
plt.savefig(os.path.join(output_directory, 'TAP_vs_RMSD_Correlation_Heatmap.png'))
plt.show()

# Perform PCA
scaler = StandardScaler()
scaled_data = scaler.fit_transform(merged_df[valid_metrics])
pca = PCA(n_components=2)  # Reduce to 2 components for visualisation
pca_results = pca.fit_transform(scaled_data)

# Plot PCA results
plt.figure(figsize=(8, 6))
plt.scatter(pca_results[:, 0], pca_results[:, 1], c=merged_df['Average RMSD'], cmap='coolwarm')
plt.colorbar(label='Average RMSD')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA of TAP Metrics with RMSD as color')
plt.savefig(os.path.join(output_directory, 'TAP_PCA_vs_RMSD.png'))
plt.show()

# Save merged data to CSV for further analysis
merged_df.to_csv(os.path.join(output_directory, 'merged_tap_rmsd_data.csv'), index=False)



# PCR loading (1)

valid_metrics = ['Total CDR Length', 'CDR Vicinity PSH Score (Kyte & Doolittle)', 
                 'CDR Vicinity PPC Score', 'CDR Vicinity PNC Score', 'SFvCSP Score']

# Standardize the data (important for PCA)
from sklearn.preprocessing import StandardScaler
X = merged_df[valid_metrics]
X_standardized = StandardScaler().fit_transform(X)

# Perform PCA
pca = PCA(n_components=2)  # PC1 and PC2
principal_components = pca.fit_transform(X_standardized)

# Create a DataFrame with the PCA results
pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

# Look at the loadings (contributions of the original features)
loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=valid_metrics)

# Plot the loadings
plt.figure(figsize=(10, 6))
plt.bar(loadings.indeximport pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os


# valid metrics for PCA 
valid_metrics = ['Total CDR Length', 'CDR Vicinity PSH Score (Kyte & Doolittle)', 
                 'CDR Vicinity PPC Score', 'CDR Vicinity PNC Score', 'SFvCSP Score']

# Standardise the data (required for PCA)
X = merged_df[valid_metrics]
X_standardized = StandardScaler().fit_transform(X)

# Perform PCA
pca = PCA(n_components=2)  # Only 2 components since you are plotting PC1 and PC2
principal_components = pca.fit_transform(X_standardized)

# Create a DataFrame with the PCA results
pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

# See loadings (contributions of the original features)
loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=valid_metrics)

# Plot the loadings
plt.figure(figsize=(10, 6))
plt.bar(loadings.index, loadings['PC1'], alpha=0.7, label='PC1')
plt.bar(loadings.index, loadings['PC2'], alpha=0.7, label='PC2', bottom=loadings['PC1'])
plt.title('Loadings of TAP Metrics on Principal Components')
plt.xlabel('TAP Metric')
plt.ylabel('Loading')
plt.legend()
plt.xticks(rotation=45)
plt.tight_layout()

# Save the plot to the output directory
plot_path = os.path.join(output_directory, 'TAP_Metrics_PCA_Loadings.png')
plt.savefig(plot_path)
plt.show()

# Display the loadings DataFrame
print(loadings)


# PCA loading (2)

# After performing PCA and have the 'pca' and 'valid_metrics' variables
loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=valid_metrics)

# Plot the loadings
plt.figure(figsize=(10, 6))
plt.bar(loadings.index, loadings['PC1'], label='PC1', alpha=0.6)
plt.bar(loadings.index, loadings['PC2'], label='PC2', alpha=0.6)
plt.xlabel('TAP Metrics')
plt.ylabel('Loading Scores')
plt.title('PCA Loadings for TAP Metrics')
plt.xticks(rotation=45, ha='right')
plt.legend()

# Save the plot to the specified directory
output_directory = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/Liability_sequence'
plot_filename = 'TAP_Metrics_PCA_Loadings.png'
plot_path = os.path.join(output_directory, plot_filename)
plt.savefig(plot_path, bbox_inches='tight')

# Show the plot
plt.show()



# PCA (3)



# Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(merged_df[valid_metrics])

# Create a DataFrame for the PCA loadings
loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=valid_metrics)

# Plot the loadings for PC1 and PC2
plt.figure(figsize=(10, 8))
loadings[['PC1', 'PC2']].plot(kind='bar', color=['#1f77b4', '#ff7f0e'], edgecolor='k')

# Annotate the bars with their respective loadings
for i in range(loadings.shape[0]):
    for j in range(loadings.shape[1]):
        plt.text(i - 0.15 + j * 0.3, loadings.iloc[i, j] + np.sign(loadings.iloc[i, j]) * 0.02,
                 f'{loadings.iloc[i, j]:.2f}', color='black', ha='center')

plt.title('PCA Loadings for TAP Metrics')
plt.ylabel('Loading Value')
plt.xlabel('TAP Metrics')
plt.xticks(rotation=45, ha='right')
plt.legend(['PC1', 'PC2'])
plt.tight_layout()

# Save the figure
output_path = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/Liability_sequence/TAP_Metrics_PCA_Loadings.png'
plt.savefig(output_path)

plt.show()
