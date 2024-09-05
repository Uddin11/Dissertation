#!/bin/bash
#SBATCH --partition=defq         # the requested queue
#SBATCH --nodes=1                # number of nodes to use
#SBATCH --ntasks=1               # number of tasks (tasks per node)
#SBATCH --cpus-per-task=64       # number of CPUs per task
#SBATCH --mem-per-cpu=8000       # memory per CPU in megabytes
#SBATCH --error=%J.err           # redirect stderr to this file
#SBATCH --output=%J.out          # redirect stdout to this file
##SBATCH --mail-user=uddinms1@Cardiff.ac.uk  # email address for notifications
##SBATCH --mail-type=end         # email on job end
##SBATCH --mail-type=fail        # email on job failure

# Load Conda (assuming conda is already in your PATH)
source /mnt/clusters/grayling/data/c22111854/ResearchProject/miniconda/etc/profile.d/conda.sh

# Activate the desired Conda environment
conda activate /mnt/clusters/grayling/data/c22111854/ResearchProject/miniconda/envs/DeepAbEnv

# Echo to verify environment activation (optional)
echo "Conda environment activated: $(which python)"

# Add the deepab directory to PYTHONPATH
export PYTHONPATH=$PYTHONPATH:/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/DeepAb

# Navigate to the directory containing predict.py
cd /mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/DeepAb/deepab

# Define the input FASTA file and the output directory
INPUT_FASTA=/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/DR1PKY15mer_L243_scFv._antibody.fasta
OUTPUT_DIR=/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/DR1PKY15mer_L243_scFv_Predictions

# Create the output directory
mkdir -p $OUTPUT_DIR

# Extract the basename of the input FASTA file
BASENAME=$(basename "$INPUT_FASTA" .fasta)

# Run the prediction
python predict.py $INPUT_FASTA --decoys 5 --renumber --pred_dir $OUTPUT_DIR --target $BASENAME

# Deactivate the Conda environment
conda deactivate

# Echo to indicate end of script (optional)
echo "Job complete"
