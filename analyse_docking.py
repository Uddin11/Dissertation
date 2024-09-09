import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define the directory containing docking outputs
docking_dir = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/L243/Docking_outputs"

# Load docking results
docked_files = [
    "fixed_DR1PKY15mer_L243_scFv._antibody.deepab.0_docked.pdb",
    "fixed_DR1PKY15mer_L243_scFv._antibody.deepab.1_docked.pdb",
    "fixed_DR1PKY15mer_L243_scFv._antibody.deepab.2_docked.pdb",
    "fixed_DR1PKY15mer_L243_scFv._antibody.deepab.3_docked.pdb",
    "fixed_DR1PKY15mer_L243_scFv._antibody.deepab.4_docked.pdb",
    "fixed_DR1PKY15mer_L243_scFv._antibody.mds_docked.pdb"
]

# Simulated data, replace this with actual calculations or loading of data
data = {
    'Complex': docked_files,
    'Total Score': [-],  # Add values
    'Interface Score (ΔG)': [.5]  # Add values
}

df = pd.DataFrame(data)

# Print out DataFrame to verify the data
print(df)

# Plotting the binding energies (Bar plot)
plt.figure(figsize=(10, 6))
sns.barplot(x='Complex', y='Interface Score (ΔG)', data=df)
plt.title('Binding Energies (ΔG) of Docked Complexes')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.join(docking_dir, 'binding_energies_comparison.png'), dpi=300)
plt.show()  # Add this to visualize the plot before saving

# Plotting the distribution of binding energies (Box plot)
plt.figure(figsize=(8, 5))
sns.boxplot(y='Interface Score (ΔG)', data=df)
plt.title('Distribution of Binding Energies (ΔG)')
plt.tight_layout()
plt.savefig(os.path.join(docking_dir, 'binding_energies_distribution.png'), dpi=300)
plt.show()

# Scatter plot of Total Score vs. Interface Score (ΔG)
plt.figure(figsize=(8, 5))
sns.scatterplot(x='Total Score', y='Interface Score (ΔG)', data=df)
plt.title('Total Score vs. Interface Score (ΔG)')
plt.tight_layout()
plt.savefig(os.path.join(docking_dir, 'score_vs_binding_energy.png'), dpi=300)
plt.show()

# Save the statistics summary to a file
stats_summary = df.describe()
stats_summary.to_csv(os.path.join(docking_dir, 'binding_energies_stats_summary.csv'))

# Display the summary table in the output
print("Summary Table:")
print(stats_summary)

print("Analysis complete. Plots and summary statistics have been saved.")

