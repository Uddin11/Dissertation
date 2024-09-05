import os
import pandas as pd

# Define CDR ranges
CDR_RANGES = {
    'H_1': (1, 10),
    'H_2': (11, 20),
    'H_3': (21, 30),
    'L_1': (1, 10),
    'L_2': (11, 20),
    'L_3': (21, 30)
}

def read_renumbered(file_path):
    """Read a renumbered file and return a dictionary with chain identifiers and sequences."""
    sequences = {}
    chain_type = None

    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue  # Skip empty lines and comments
                if line.startswith('H ') or line.startswith('L '):
                    parts = line.split()
                    if len(parts) > 1:
                        chain_id = parts[0]  # Chain identifier (e.g., 'H', 'L')
                        position = int(parts[1])  # Position
                        residue = parts[2]  # Amino acid residue

                        if chain_id not in sequences:
                            sequences[chain_id] = []
                        
                        sequences[chain_id].append(residue)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")

    # Join sequences into single strings
    for chain_id in sequences:
        sequences[chain_id] = ''.join(sequences[chain_id])
    
    return sequences

def extract_cdr_region(sequence, start, end):
    """Extract the CDR region from the sequence."""
    if sequence:  # Check if sequence is not empty
        return sequence[start-1:end]  # Convert to 0-based index
    else:
        return ''

def process_directory(directory, cdr_ranges):
    """Process all renumbered files in a directory and extract CDR regions."""
    cdr_data = []
    for filename in os.listdir(directory):
        if filename.endswith('.renumbered'):
            file_path = os.path.join(directory, filename)
            if os.path.exists(file_path):
                print(f"Reading file: {file_path}")
                renumbered_sequences = read_renumbered(file_path)
                for chain_id, sequence in renumbered_sequences.items():
                    # Determine if chain_id is H or L and extract corresponding CDR regions
                    chain_type = 'H' if chain_id.startswith('H') else 'L'
                    for cdr_name, (start, end) in cdr_ranges.items():
                        if cdr_name.startswith(chain_type):
                            cdr_sequence = extract_cdr_region(sequence, start, end)
                            cdr_data.append({
                                'File': filename,
                                'Sequence ID': chain_id,
                                'CDR Region': cdr_name,
                                'Sequence': cdr_sequence
                            })
            else:
                print(f"File not found: {file_path}")
    return cdr_data

def main():
    # Define input directories
    bound_dir = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/bound_renumbered'
    unbound_dir = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/unbound_renumbered'
    
    # Process both directories
    bound_cdr_data = process_directory(bound_dir, CDR_RANGES)
    unbound_cdr_data = process_directory(unbound_dir, CDR_RANGES)
    
    # Combine data from both directories
    all_cdr_data = bound_cdr_data + unbound_cdr_data
    
    # Save extracted CDR regions to CSV
    cdr_df = pd.DataFrame(all_cdr_data)
    cdr_csv_path = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/extracted_cdr_regions.csv'
    cdr_df.to_csv(cdr_csv_path, index=False)
    print(f"Extracted CDR regions saved to {cdr_csv_path}")

if __name__ == '__main__':
    main()



#####


import pandas as pd

# Path to the filtered CSV file
filtered_csv_file_path = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/filtered_sequences.csv'

# Load the filtered DataFrame
filtered_df = pd.read_csv(filtered_csv_file_path)

# List of required CDR regions
required_regions = ['H_1', 'H_2', 'H_3', 'L_1', 'L_2', 'L_3']

# Function to check completeness
def check_completeness(df):
    missing_data = []
    for file_id, group in df.groupby('File'):
        missing_regions = set(required_regions) - set(group['CDR Region'])
        if missing_regions:
            missing_data.append((file_id, missing_regions))
    return missing_data

# Check for missing CDR regions
missing_data = check_completeness(filtered_df)

# Display missing data
if missing_data:
    print("Antibodies with missing CDR regions:")
    for file_id, missing_regions in missing_data:
        print(f"{file_id} is missing: {', '.join(missing_regions)}")
else:
    print("All antibodies of interest have all required CDR regions.")


#######################################################


# 'Padding' the distances so they are the same length. Hamming metric requires sequences to be the same lenth 

import pandas as pd
from scipy.spatial.distance import hamming

# Load your filtered DataFrame
def load_filtered_df(file_path):
    return pd.read_csv(file_path)

# Function to pad or trim sequences to the same length
def pad_or_trim_sequence(sequence, length):
    if len(sequence) < length:
        return sequence + '-' * (length - len(sequence))
    elif len(sequence) > length:
        return sequence[:length]
    else:
        return sequence

# Function to get CDR sequences for a given antibody ID
def get_cdr_sequence(df, file_id):
    df_subset = df[df['File'].str.startswith(file_id)]
    if df_subset.empty:
        return None
    return ''.join(df_subset.sort_values('CDR Region')['Sequence'])

# Calculate the Hamming distance between two sequences
def calculate_hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length for Hamming distance.")
    return hamming(list(seq1), list(seq2))

# Main function to perform the comparison
def find_top_matches(filtered_df, antibodies_of_interest):
    # Extract CDR sequences and determine maximum length for padding
    cdr_sequences = {file_id: get_cdr_sequence(filtered_df, file_id) for file_id in filtered_df['File'].unique()}
    cdr_sequences = {k: v for k, v in cdr_sequences.items() if v is not None}  # Filter out None values
    
    if not cdr_sequences:
        raise ValueError("No valid CDR sequences found. Please check your data or filtering criteria.")
    
    max_length = max(len(seq) for seq in cdr_sequences.values())
    
    # Pad or trim sequences
    cdr_sequences_padded = {file_id: pad_or_trim_sequence(seq, max_length) for file_id, seq in cdr_sequences.items()}

    results = {}
    for interest_id in antibodies_of_interest:
        interest_seq = pad_or_trim_sequence(get_cdr_sequence(filtered_df, interest_id), max_length)
        
        if not interest_seq:
            print(f"Warning: No sequence found for {interest_id}. Skipping.")
            continue
        
        print(f"Comparing {interest_id} with sequences:")
        distances = []
        for other_id, other_seq in cdr_sequences_padded.items():
            if other_id != interest_id:
                # Calculate Hamming distance
                distance = calculate_hamming_distance(interest_seq, other_seq)
                distances.append((other_id, distance))
        
        # Find the top 10 closest matches
        distances.sort(key=lambda x: x[1])
        top_matches = distances[:10]
        results[interest_id] = top_matches
        print(f"Top matches for {interest_id}: {top_matches}")

    return results

# Display results
def display_results(results):
    for interest_id, matches in results.items():
        print(f"\nTop 10 matches for {interest_id}:")
        for match_id, dist in matches:
            print(f"ID: {match_id}, Hamming Distance: {dist:.4f}")

if __name__ == "__main__":
    # Path to the filtered CSV file
    filtered_csv_file_path = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/extracted_cdr_regions.csv'
    
    # Load filtered DataFrame
    filtered_df = load_filtered_df(filtered_csv_file_path)
    
    # List of antibodies of interest
    antibodies_of_interest = [
        '1cfs', '1cfv', '1clz', '1cu4', '1d6v', '1dbj', '1dbm',
        '1a2y', '1a3r', '1a6v', '1acy', '1afv', '1ai1', '1ar1',
        '1baf', '1bgx', '1bog', '1bvk', '1c12', '1c5c', '1ce1',
        '1cfn', '1cft', '1cly', '1ct8', '1cz8', '1dbb', '1dbk',
        '1a14', '1a3l', '1a4k', '1a6w', '1adq', '1ahw', '1aj7',
        '1axs', '1bfv', '1bj1', '1bql', '1c08', '1c1e', '1cbv',
        '1cf8'
    ]

    # Find and display top matches
    results = find_top_matches(filtered_df, antibodies_of_interest)
    display_results(results)



# Hamming Distamce and generates visuals 

import pandas as pd
from scipy.spatial.distance import hamming
import matplotlib.pyplot as plt
import os

# Load your filtered DataFrame
def load_filtered_df(file_path):
    return pd.read_csv(file_path)

# Function to pad or trim sequences to the same length
def pad_or_trim_sequence(sequence, length):
    if len(sequence) < length:
        return sequence + '-' * (length - len(sequence))
    elif len(sequence) > length:
        return sequence[:length]
    else:
        return sequence

# Function to get CDR sequences for a given antibody ID
def get_cdr_sequence(df, file_id):
    df_subset = df[df['File'].str.startswith(file_id)]
    if df_subset.empty:
        return None
    return ''.join(df_subset.sort_values('CDR Region')['Sequence'])

# Calculate the Hamming distance between two sequences
def calculate_hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length for Hamming distance.")
    return hamming(list(seq1), list(seq2))

# Main function to perform the comparison
def find_top_matches(filtered_df, antibodies_of_interest):
    # Extract CDR sequences and determine maximum length for padding
    cdr_sequences = {file_id: get_cdr_sequence(filtered_df, file_id) for file_id in filtered_df['File'].unique()}
    cdr_sequences = {k: v for k, v in cdr_sequences.items() if v is not None}  # Filter out None values
    
    if not cdr_sequences:
        raise ValueError("No valid CDR sequences found. Please check your data or filtering criteria.")
    
    max_length = max(len(seq) for seq in cdr_sequences.values())
    
    # Pad or trim sequences
    cdr_sequences_padded = {file_id: pad_or_trim_sequence(seq, max_length) for file_id, seq in cdr_sequences.items()}

    results = {}
    for interest_id in antibodies_of_interest:
        interest_seq = pad_or_trim_sequence(get_cdr_sequence(filtered_df, interest_id), max_length)
        
        if not interest_seq:
            print(f"Warning: No sequence found for {interest_id}. Skipping.")
            continue
        
        print(f"Comparing {interest_id} with sequences:")
        distances = []
        for other_id, other_seq in cdr_sequences_padded.items():
            if other_id != interest_id:
                # Calculate Hamming distance
                distance = calculate_hamming_distance(interest_seq, other_seq)
                distances.append((other_id, distance))
        
        # Find the top 10 closest matches
        distances.sort(key=lambda x: x[1])
        top_matches = distances[:10]
        results[interest_id] = top_matches
        print(f"Top matches for {interest_id}: {top_matches}")

    return results

import os
import matplotlib.pyplot as plt

# Display and save the results as images
def plot_and_save_results(results, output_dir='/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures'):
    # Create the output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Iterate over the results to create plots
    for interest_id, matches in results.items():
        ids = [match[0] for match in matches]
        distances = [match[1] for match in matches]
        
        # Create a horizontal bar plot
        plt.figure(figsize=(10, 6))
        plt.barh(ids, distances, color='blue')
        plt.xlabel('Hamming Distance')
        plt.title(f'Top 10 Matches for {interest_id}')
        plt.gca().invert_yaxis()  # Invert y-axis to have the lowest distance at the top
        plt.grid(axis='x', linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        # Save the figure as a PNG file
        plt.savefig(os.path.join(output_dir, f'{interest_id}_top_matches.png'))
        plt.close()

if __name__ == "__main__":
    # Path to the filtered CSV file
    filtered_csv_file_path = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/extracted_cdr_regions.csv'
    
    # Load filtered DataFrame
    filtered_df = load_filtered_df(filtered_csv_file_path)
    
    # List of antibodies of interest
    antibodies_of_interest = [
        '1cfs', '1cfv', '1clz', '1cu4', '1d6v', '1dbj', '1dbm',
        '1a2y', '1a3r', '1a6v', '1acy', '1afv', '1ai1', '1ar1',
        '1baf', '1bgx', '1bog', '1bvk', '1c12', '1c5c', '1ce1',
        '1cfn', '1cft', '1cly', '1ct8', '1cz8', '1dbb', '1dbk',
        '1a14', '1a3l', '1a4k', '1a6w', '1adq', '1ahw', '1aj7',
        '1axs', '1bfv', '1bj1', '1bql', '1c08', '1c1e', '1cbv',
        '1cf8'
    ]

    # Find top matches
    results = find_top_matches(filtered_df, antibodies_of_interest)
    
    # Plot and save results
    plot_and_save_results(results)


# puts visuals into one image

import os
import matplotlib.pyplot as plt
import math

# Display and save all the results as a single image
def plot_and_save_all_results(results, output_file='/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/all_top_matches.png'):
    num_plots = len(results)
    cols = 4  # Number of columns in the grid
    rows = math.ceil(num_plots / cols)  # Calculate number of rows required

    fig, axes = plt.subplots(rows, cols, figsize=(20, 5 * rows))
    axes = axes.flatten()  # Flatten the 2D array of axes into a 1D array for easier indexing

    for ax, (interest_id, matches) in zip(axes, results.items()):
        ids = [match[0] for match in matches]
        distances = [match[1] for match in matches]
        
        ax.barh(ids, distances, color='blue')
        ax.set_xlabel('Hamming Distance')
        ax.set_title(f'Top 10 Matches for {interest_id}')
        ax.invert_yaxis()
        ax.grid(axis='x', linestyle='--', alpha=0.7)
    
    # Hide any unused subplots
    for i in range(num_plots, len(axes)):
        fig.delaxes(axes[i])

    plt.tight_layout()
    
    # Save the figure as a single PNG file
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    # Path to the filtered CSV file
    filtered_csv_file_path = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/extracted_cdr_regions.csv'
    
    # Load filtered DataFrame
    filtered_df = load_filtered_df(filtered_csv_file_path)
    
    # List of antibodies of interest
    antibodies_of_interest = [
        '1cfs', '1cfv', '1clz', '1cu4', '1d6v', '1dbj', '1dbm',
        '1a2y', '1a3r', '1a6v', '1acy', '1afv', '1ai1', '1ar1',
        '1baf', '1bgx', '1bog', '1bvk', '1c12', '1c5c', '1ce1',
        '1cfn', '1cft', '1cly', '1ct8', '1cz8', '1dbb', '1dbk',
        '1a14', '1a3l', '1a4k', '1a6w', '1adq', '1ahw', '1aj7',
        '1axs', '1bfv', '1bj1', '1bql', '1c08', '1c1e', '1cbv',
        '1cf8'
    ]

    # Find top matches
    results = find_top_matches(filtered_df, antibodies_of_interest)
    
    # Plot and save all results in one image
    plot_and_save_all_results(results)


# combine the plots 


from PIL import Image
import os
import math

def combine_images(image_files, output_file, images_per_row=4):
    # Load all images
    images = [Image.open(img_file) for img_file in image_files]
    
    # Determine grid size
    num_images = len(images)
    num_rows = math.ceil(num_images / images_per_row)
    
    # Determine size of each image
    img_width, img_height = images[0].size
    
    # Create a new blank image with the appropriate size
    combined_width = images_per_row * img_width
    combined_height = num_rows * img_height
    combined_image = Image.new('RGB', (combined_width, combined_height), 'white')
    
    # Paste each image into the combined image
    for index, img in enumerate(images):
        row = index // images_per_row
        col = index % images_per_row
        x = col * img_width
        y = row * img_height
        combined_image.paste(img, (x, y))
    
    # Save the combined image
    combined_image.save(output_file)

if __name__ == "__main__":
    # Directory containing the individual PNG files
    input_dir = '/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures'
    
    # Output file for the combined image
    output_file = os.path.join(input_dir, 'combined_top_matches.png')
    
    # List of individual PNG files to combine
    image_files = sorted([os.path.join(input_dir, file) for file in os.listdir(input_dir) if file.endswith('_top_matches.png')])
    
    # Combine the images into one
    combine_images(image_files, output_file)
