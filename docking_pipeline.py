import pyrosetta
from pyrosetta import init, pose_from_file, Pose
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

# Initialize PyRosetta
init()

# Set up logging
logging.basicConfig(filename='docking_pipeline.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Define paths
protein_pdb = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/L243/fixed_8pje.pdb"
antibody_pdbs = [
    "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/DR1PKY15mer_L243_scFv_Predictions/decoys_pdb_fixer_output/fixed_DR1PKY15mer_L243_scFv._antibody.deepab.0.pdb",
    "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/DR1PKY15mer_L243_scFv_Predictions/decoys_pdb_fixer_output/fixed_DR1PKY15mer_L243_scFv._antibody.deepab.1.pdb",
    "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/DR1PKY15mer_L243_scFv_Predictions/decoys_pdb_fixer_output/fixed_DR1PKY15mer_L243_scFv._antibody.deepab.2.pdb",
    "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/DR1PKY15mer_L243_scFv_Predictions/decoys_pdb_fixer_output/fixed_DR1PKY15mer_L243_scFv._antibody.deepab.3.pdb",
    "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/DR1PKY15mer_L243_scFv_Predictions/decoys_pdb_fixer_output/fixed_DR1PKY15mer_L243_scFv._antibody.deepab.4.pdb",
    "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/DR1PKY15mer_L243_scFv_Predictions/decoys_pdb_fixer_output/fixed_DR1PKY15mer_L243_scFv._antibody.mds.pdb"
]

output_dir = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/L243/Docking_outputs"

def process_antibody(antibody_pdb):
    try:
        logging.info(f"Processing {antibody_pdb}...")
        
        # Load the antibody structure
        antibody_pose = pose_from_file(antibody_pdb)

        # Create a combined pose
        combined_pose = Pose()
        combined_pose.assign(protein_pose)

        # Append the antibody pose by jump
        combined_pose.append_pose_by_jump(antibody_pose, 1)

        # Define output path
        output_pdb = os.path.join(output_dir, os.path.basename(antibody_pdb).replace(".pdb", "_docked.pdb"))

        # Save the combined pose
        combined_pose.dump_pdb(output_pdb)
        logging.info(f"Docked structure saved to {output_pdb}")
        return output_pdb
    except Exception as e:
        logging.error(f"Error processing {antibody_pdb}: {e}")
        return None

def evaluate_structure(docked_pdb):
    try:
        logging.info(f"Evaluating {docked_pdb}...")
        antibody_pose = pose_from_file(docked_pdb)
        combined_pose = antibody_pose.clone()
        combined_pose.append_pose_by_jump(antigen_pose, 1)

        # Compute the binding energy
        score = pyrosetta.get_fa_scorefxn()(combined_pose)
        interface_analyzer = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover("A_B")
        interface_analyzer.apply(combined_pose)
        interface_score = interface_analyzer.get_interface_dG()
        logging.info(f"{docked_pdb}: Total Score = {score}, Interface Score (ΔG) = {interface_score}")
    except Exception as e:
        logging.error(f"Error evaluating {docked_pdb}: {e}")

# Load the protein structure
try:
    protein_pose = pose_from_file(protein_pdb)
except Exception as e:
    logging.error(f"Failed to load protein structure: {e}")
    raise

# Parallel docking process
with ProcessPoolExecutor() as executor:
    future_to_pdb = {executor.submit(process_antibody, pdb): pdb for pdb in antibody_pdbs}
    docked_pdbs = []
    for future in as_completed(future_to_pdb):
        result = future.result()
        if result:
            docked_pdbs.append(result)

logging.info("Docking process complete.")

# Parallel binding energy evaluation process
try:
    antigen_pose = pose_from_file(protein_pdb)
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(evaluate_structure, pdb) for pdb in docked_pdbs]
        for future in as_completed(futures):
            future.result()  # This will raise any exception caught during the execution
    logging.info("Binding energy calculations complete.")
except Exception as e:
    logging.error(f"Failed to perform binding energy calculations: {e}")



# plots between original predicted structure and decoys 


# Adjust the y-axis to match the scale of the data
plt.ylim(-60, 0)

plt.legend()

# Save and show the plot
output_path = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/L243/Docking_outputs/visual_inspection_plot_adjusted.png"
plt.savefig(output_path)
plt.show()


###



# Test to extract/calculate energies 

import pyrosetta
from pyrosetta import pose_from_file
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import logging

# Initialize PyRosetta
pyrosetta.init()

# Define paths
output_dir = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/L243/Docking_outputs"
summary_csv_path = os.path.join(output_dir, "binding_energies_stats_summary.csv")
score_plot_path = os.path.join(output_dir, "binding_energies_distribution.png")
score_vs_binding_plot_path = os.path.join(output_dir, "score_vs_binding_energy.png")

# Ensure logging is configured
logging.basicConfig(filename='docking_pipeline.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def collect_results():
    results = []
    pdb_files = [
        "fixed_DR1PKY15mer_L243_scFv._antibody.deepab.0_docked.pdb",
        "fixed_DR1PKY15mer_L243_scFv._antibody.deepab.1_docked.pdb"  ## add the rest of the files here 
    ]
    for pdb in pdb_files:
        pdb_path = os.path.join(output_dir, pdb)
        try:
            # Load the docked pose
            pose = pose_from_file(pdb_path)
            # Compute the total score
            score = pyrosetta.get_fa_scorefxn()(pose)
            # Compute the interface score (ΔG)
            interface_analyzer = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover("A_B")
            interface_analyzer.apply(pose)
            interface_score = interface_analyzer.get_interface_dG()
            results.append({
                "PDB File": os.path.basename(pdb),
                "Total Score": score,
                "Interface Score (ΔG)": interface_score
            })
        except Exception as e:
            logging.error(f"Error collecting results for {pdb_path}: {e}")
    
    # Save results to CSV
    df_results = pd.DataFrame(results)
    df_results.to_csv(summary_csv_path, index=False)
    logging.info(f"Docking results saved to {summary_csv_path}")
    
    return df_results

def generate_plots(df_results):
    # Distribution of total scores
    plt.figure(figsize=(10, 6))
    sns.histplot(df_results['Total Score'], kde=True, color='skyblue')
    plt.xlabel('Total Score')
    plt.ylabel('Frequency')
    plt.title('Distribution of Docking Scores')
    plt.tight_layout()
    plt.savefig(score_plot_path)
    plt.close()
    
    # Distribution of interface scores
    plt.figure(figsize=(10, 6))
    sns.histplot(df_results['Interface Score (ΔG)'], kde=True, color='salmon')
    plt.xlabel('Interface Score (ΔG)')
    plt.ylabel('Frequency')
    plt.title('Distribution of Interface Scores')
    plt.tight_layout()
    plt.savefig(score_vs_binding_plot_path)
    plt.close()
    
    logging.info(f"Plots saved to {score_plot_path} and {score_vs_binding_plot_path}")

# Main execution
if __name__ == "__main__":
    df_results = collect_results()
    generate_plots(df_results)
    logging.info("Analysis and visualization complete.")


