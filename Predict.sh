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

set -e  # Exit immediately if a command exits with a non-zero status

# Load Conda (assuming conda is already in your PATH)
source /mnt/clusters/grayling/data/c22111854/ResearchProject/miniconda/etc/profile.d/conda.sh

# Activate the desired Conda environment
conda activate /mnt/clusters/grayling/data/c22111854/ResearchProject/miniconda/envs/DeepAbEnv

# Verify conda environment activation
if [[ $? -ne 0 ]]; then
  echo "Failed to activate conda environment" >&2
  exit 1
fi

# Ensure the correct Python environment is being used
export PATH="/mnt/clusters/grayling/data/c22111854/ResearchProject/miniconda/envs/DeepAbEnv/bin:$PATH"
export PYTHONPATH="/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/DeepAb:$PYTHONPATH"
export PYTORCH_JIT=0

# Check the Python executable path
echo "Using Python from: $(which python)"

# Check PYTHONPATH
echo "PYTHONPATH: $PYTHONPATH"

# Check installed packages
python -m pip list | grep tqdm

# Check if tqdm is installed
python -c "import tqdm" 2>/dev/null || { echo 'tqdm is not installed' >&2; exit 1; }

# Define input directories
input_dirs=(
  "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/bound"
  "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/unbound"
)

# Define output base directory
output_base_dir="/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/SabDab_DeepAb_Prediction3"

# Create output base directory if it doesn't exist
mkdir -p $output_base_dir

# Loop over each input directory
for input_dir in "${input_dirs[@]}"; do
    # Determine the subdirectory name based on whether it's bound or unbound
    if [[ $input_dir == *"bound"* ]]; then
        sub_dir="bound"
    else
        sub_dir="unbound"
    fi

    # Create the subdirectory within the output base directory
    output_dir="$output_base_dir/$sub_dir"
    mkdir -p $output_dir

    # Loop over each FASTA file in the input directory and run prediction
    for fasta_file in $input_dir/*.fasta; do
        if [[ -f $fasta_file ]]; then
            # Extract the basename without the extension
            base_name=$(basename $fasta_file .fasta)
            
            # Create a subdirectory for the output
            fasta_output_dir="$output_dir/$base_name"
            mkdir -p $fasta_output_dir
            
            # Run prediction for the current fasta file
            python /mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/DeepAb/deepab/predict.py \
                $fasta_file \
                --pred_dir $fasta_output_dir \
                --decoys 5 \
                --renumber \
                --use_gpu
        else
            echo "No FASTA files found in the input directory $input_dir" >&2
            exit 1
        fi
    done
done

# Deactivate the Conda environment
conda deactivate

# Echo to indicate end of script (optional)
echo "Job complete"
