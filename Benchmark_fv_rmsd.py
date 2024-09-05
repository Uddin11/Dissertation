import Bio.PDB
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
from sklearn.manifold import TSNE
from PIL import Image

def calculate_rmsd(native_structure, decoy_structure):
    super_imposer = Bio.PDB.Superimposer()
    native_fv = get_fv_region(native_structure)
    pred_chains = list(decoy_structure.get_chains())
    ref_atoms = []
    alt_atoms = []
    for pred_chain, native_chain in zip(pred_chains, native_fv):
        for pred_res, native_res in zip(pred_chain, native_chain):
            if pred_res.id[1] == native_res.id[1]:  # Match residues by their sequence number
                if 'CA' in native_res and 'CA' in pred_res:
                    ref_atoms.append(native_res['CA'])
                    alt_atoms.append(pred_res['CA'])
    if len(ref_atoms) == 0 or len(alt_atoms) == 0:
        raise ValueError("No matching residues found for RMSD calculation.")
    super_imposer.set_atoms(ref_atoms, alt_atoms)
    super_imposer.apply(decoy_structure.get_atoms())
    return super_imposer.rms

def get_fv_region(structure):
    fv_chains = [chain for chain in structure.get_chains() if chain.id in ['H', 'L']]
    return fv_chains

def analyze_rmsd(rmsd_values, title):
    rmsd_array = np.array(rmsd_values)
    mean_rmsd = np.mean(rmsd_array)
    median_rmsd = np.median(rmsd_array)
    std_dev_rmsd = np.std(rmsd_array)
    
    print(f"\n{title} - Mean RMSD: {mean_rmsd:.3f}")
    print(f"{title} - Median RMSD: {median_rmsd:.3f}")
    print(f"{title} - Standard Deviation of RMSD: {std_dev_rmsd:.3f}")

    # Create a figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Bar chart
    axes[0].bar(range(len(rmsd_array)), rmsd_array, color='blue')
    axes[0].set_xlabel('Decoy Index')
    axes[0].set_ylabel('RMSD')
    axes[0].set_title(f'{title} - RMSD Bar Chart')

    # Heatmap
    rmsd_matrix = np.zeros((len(rmsd_array), len(rmsd_array)))
    for i in range(len(rmsd_array)):
        for j in range(len(rmsd_array)):
            rmsd_matrix[i, j] = np.abs(rmsd_array[i] - rmsd_array[j])
    sns.heatmap(rmsd_matrix, ax=axes[1], cmap="viridis", annot=True, xticklabels=range(len(rmsd_array)), yticklabels=range(len(rmsd_array)))
    axes[1].set_xlabel('Decoy Index')
    axes[1].set_ylabel('Decoy Index')
    axes[1].set_title(f'{title} - RMSD Heatmap')

    # Density (KDE) plot - Check for variance
    if std_dev_rmsd > 0:
        sns.kdeplot(rmsd_array, ax=axes[2], fill=True)
        axes[2].set_xlabel('RMSD')
        axes[2].set_ylabel('Density')
        axes[2].set_title(f'{title} - RMSD KDE Plot')
    else:
        axes[2].text(0.5, 0.5, 'No variance in RMSD values', ha='center', va='center', fontsize=12)
        axes[2].set_title(f'{title} - No Variance')
        axes[2].axis('off')

    # Adjust layout and save all plots in one PNG file in the specified directory
    plt.tight_layout()
    plt.savefig(f'/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/{title}_RMSD_analysis.png')
    plt.show()

def process_directory(base_native_path, base_decoy_path, directory):
    native_pdb_path = os.path.join(base_native_path, f"{directory}.pdb")
    
    if not os.path.isfile(native_pdb_path):
        print(f"Native file not found for {directory}: {native_pdb_path}")
        return
    
    native_structure = pdb_parser.get_structure("native", native_pdb_path)
    rmsd_results = []

    for subdirectory in os.listdir(os.path.join(base_decoy_path, directory)):
        subdirectory_path = os.path.join(base_decoy_path, directory, subdirectory)
        if not os.path.isdir(subdirectory_path):
            continue
        
        decoy_paths = [os.path.join(subdirectory_path, f) for f in os.listdir(subdirectory_path) if f.endswith('.pdb')]

        for decoy_path in decoy_paths:
            try:
                decoy_structure = pdb_parser.get_structure("decoy", decoy_path)
                rmsd = calculate_rmsd(native_structure, decoy_structure)
                rmsd_results.append((os.path.basename(decoy_path), rmsd))
                print(f"{directory}/{subdirectory} - {os.path.basename(decoy_path)} RMSD: {rmsd:.3f}")
            except ValueError as ve:
                print(f"ValueError processing {decoy_path}: {ve}")
            except Exception as e:
                print(f"Error processing {decoy_path}: {e}")

    if rmsd_results:
        rmsd_values = [result[1] for result in rmsd_results]
        # Save results to CSV
        csv_filename = f"/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/{directory}_rmsd_results.csv"
        df = pd.DataFrame(rmsd_results, columns=["Decoy Name", "RMSD Value"])
        df.to_csv(csv_filename, index=False)
        print(f"RMSD results saved to {csv_filename}")
        # Analyze RMSD
        analyze_rmsd(rmsd_values, f"{directory} RMSD Analysis")
        # Generate t-SNE plot for each directory
        generate_tsne_plot(df, directory)
        # Return RMSD values for further analysis
        return df
    return None

def generate_tsne_plot(df, title):
    # We assume the DataFrame df contains RMSD values which will be used for t-SNE
    features = df[['RMSD Value']].values  # Converting to numpy array
    tsne = TSNE(n_components=2, perplexity=30, n_iter=1000, random_state=42)
    tsne_results = tsne.fit_transform(features)

    plt.figure(figsize=(10, 8))
    plt.scatter(tsne_results[:, 0], tsne_results[:, 1], c='blue', cmap='viridis')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title(f'{title} - t-SNE Visualization')
    plt.savefig(f'/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/{title}_tsne_analysis.png')
    plt.show()

def compare_rmsd(df, directory):
    original_rmsd = df[df["Decoy Name"] == "pred.deepab.0.pdb"]["RMSD Value"].values[0]
    lower_rmsd_count = sum(df[df["Decoy Name"] != "pred.deepab.0.pdb"]["RMSD Value"] < original_rmsd)
    return directory, lower_rmsd_count

def main():
    global pdb_parser
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)

    base_native_path = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/chothia/"
    base_decoy_path = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/SabDab_DeepAb_Prediction3/bound/"

    directories = [d for d in os.listdir(base_decoy_path) if os.path.isdir(os.path.join(base_decoy_path, d))]

    all_results = []

    for directory in directories:
        df = process_directory(base_native_path, base_decoy_path, directory)
        if df is not None:
            comparison_result = compare_rmsd(df, directory)
            all_results.append(comparison_result)

    if all_results:
        # Plot the comparison results
        result_df = pd.DataFrame(all_results, columns=["Antibody", "Lower RMSD Count"])
        result_df = result_df.sort_values(by="Lower RMSD Count", ascending=False)

        plt.figure(figsize=(10, 8))
        sns.barplot(x="Lower RMSD Count", y="Antibody", data=result_df, palette="viridis")
        plt.xlabel("Number of Lower RMSD Predictions (1-5)")
        plt.ylabel("Antibody")
        plt.title("Count of Predictions with Lower RMSD than Original Structure")
        plt.tight_layout()
        plt.savefig('/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/lower_rmsd_better_visualization.png')
        plt.show()

if __name__ == "__main__":
    main()

# Collate all above PNGs

# Directory where individual images are stored
image_directory = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/'
image_files = [f for f in os.listdir(image_directory) if f.endswith('_RMSD_analysis.png')]

# Open images and store in a list
images = [Image.open(os.path.join(image_directory, img)) for img in image_files]

# Calculate total width and max height for the final image
widths, heights = zip(*(img.size for img in images))
total_width = sum(widths)
max_height = max(heights)

# Create a new blank image with calculated width and height
collated_image = Image.new('RGB', (total_width, max_height))

# Paste individual images side by side
x_offset = 0
for img in images:
    collated_image.paste(img, (x_offset, 0))
    x_offset += img.width

# Save the final collated image in the specified directory
collated_image.save(os.path.join(image_directory, 'collated_rmsd_analysis.png'))

# Summary of RMSDs

csv_directory = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/'
csv_files = [f for f in os.listdir(csv_directory) if f.endswith('_rmsd_results.csv')]

average_rmsds = []
antibody_names = []

for csv_file in csv_files:
    df = pd.read_csv(os.path.join(csv_directory, csv_file))
    average_rmsd = df['RMSD Value'].mean()
    antibody_name = csv_file.split('_')[0]
    average_rmsds.append(average_rmsd)
    antibody_names.append(antibody_name)

plt.figure(figsize=(12, 8))
plt.bar(antibody_names, average_rmsds, color='skyblue')
plt.xlabel('Antibody')
plt.ylabel('Average RMSD')
plt.title('Average RMSD for Each Antibody')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(csv_directory, 'average_rmsd_overview.png'))
plt.show()
