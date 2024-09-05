#!/bin/bash
#SBATCH --job-name=sc_pipeline
#SBATCH --output=sc_pipeline_%j.out
#SBATCH --error=sc_pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000

echo "Some Usable Environment Variables:"
echo "================================="
echo "hostname=$(hostname)"
echo "\$SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "\$SLURM_NTASKS=${SLURM_NTASKS}"
echo "\$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}"
echo "\$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "\$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}"
echo "\$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}"

# Source the CCP4 environment setup script
echo "Sourcing CCP4 setup script..."
source /mnt/clusters/grayling/data/c22111854/ResearchProject/ccp4-9/bin/ccp4.setup-sh
echo "CCP4 setup script sourced."

# Source Conda environment activation script
echo "Sourcing Conda setup script..."
source /mnt/clusters/grayling/data/c22111854/ResearchProject/miniconda/etc/profile.d/conda.sh
echo "Conda setup script sourced."

# Activate the desired Conda environment
echo "Activating Conda environment..."
conda activate /mnt/clusters/grayling/data/c22111854/ResearchProject/miniconda/envs/DeepAbEnv
echo "Conda environment activated."

# Navigate to the directory where your Python script is located
cd /mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/Scripts

# Run your Python script
echo "Running Python script..."
python sc.py
echo "Python script finished."

# Additional debug information
echo "Listing current directory contents..."
ls -l
echo "Current directory contents listed."
