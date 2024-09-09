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

# Run the Python docking script with debugging enabled
python /mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/Scripts/analyse_docking.py

# Deactivate the Conda environment
conda deactivate

echo "Job complete"
