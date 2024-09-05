import numpy as np
import pandas as pd
import glob
import multiprocessing as mp
from joblib import Parallel, delayed
import os

# Define a class to handle immunoglobulin (antibody) data
class Immunoglobulin:
    total = 0  # Class variable to keep track of the number of instances

    def __init__(self, path, im_chains, ag_chains, outfile):
        Immunoglobulin.total += 1  # Increment the total count of instances
        self.path = path
        self.name = path.split("/")[-1].replace(".pdb", "")  # Extract file name without extension
        self.outfile = outfile

        self.immuno_chains = self.get_chains(im_chains)  # Get immunoglobulin chains
        self.antigen_chains = self.get_chains(ag_chains)  # Get antigen chains

        # Define the CDR ranges
        self.cdr_ranges = {"CDR1": range(27, 39),
                           "CDR2": range(56, 66),
                           "CDR3": range(105, 118)}

        # Define the antigen ranges
        self.antigen_ranges = {"antigen": range(0, 10000)}

        # Check if the file has too many atoms
        self.too_many_atoms = self.count_atoms()

        if self.too_many_atoms:
            print(f"Skipping {self.name}: too many atoms")
        else:
            self.immuno_pdb = self.retrieve_pdb_lines(self.immuno_chains, self.cdr_ranges)  # Retrieve PDB lines for immunoglobulin chains
            self.antigen_pdb = self.retrieve_pdb_lines(self.antigen_chains, self.antigen_ranges)  # Retrieve PDB lines for antigen chains
            self.distances = None  # Set distances as None

    def count_atoms(self):
        with open(self.path) as f:
            for line in f:
                if line.startswith("ATOM") and int(line[6:11]) >= 100000:  # Check if any atom number is >= 100000
                    return True
        return False

    def get_chains(self, chains):
        data = {}
        with open(self.path) as f:  # Open the PDB file
            for line in f:
                if line.startswith("REMARK") and "=" in line:  # Look for lines with REMARK and "="
                    line = line.split()  # Splits the line into a list of substrings
                    for chain in chains:  # Iterates over each chain specified in the `chains` dictionary
                        chainstr = chains[chain]  # Retrieves the chain identifier string from the `chains` dictionary
                        index = list(filter(lambda x: x.startswith(chainstr), line))  # Filters the split line to find substrings starting with `chainstr`
                        if index:  # If `index` is not empty
                            residues = index[0].split("=")[-1]  # Extracts the substring after the last "=" sign
                            if residues != "NONE":  # If residues are not "NONE"
                                residues = residues.split(";")  # Splits residues by ";" into a list
                                data[chain] = residues  # Stores the list of residues in the `data` dictionary
                                
        print(f"Chains found in {self.name}: {data}")  # Prints the chains found
        return data

    def retrieve_residues(self, c, res):
        lines = []
        amino_acids = ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'GLY', 'HIS', 'LEU',
                       'ARG', 'TRP', 'ILE', 'PRO', 'THR', 'PHE', 'ASN', 'ALA',
                       'VAL', 'GLU', 'TYR', 'MET']

        with open(self.path) as f:
            for line in f:
                if line.startswith("ATOM"):  # Look for the line starting with "ATOM"
                    chain = line[21].strip() == c  # Check if the chain matches
                    residues = int(line[22:26].strip()) in list(res)  # Check if the residue is in the list
                    amino_acid = line[17:20].strip() in amino_acids  # Check if it's a recognized amino acid
                    if chain and residues and amino_acid:
                        lines.append(line)  # If all checks pass, add the line to the list
        print(f"Residues found in chain {c} of {self.name}: {len(lines)}")
        return lines  # Return the list of lines

    def retrieve_pdb_lines(self, chains_of_interest, range_of_interest):
        pdb_chunks = {}

        if chains_of_interest:
            for chain in chains_of_interest:
                for cdr in range_of_interest:
                    subchains = chains_of_interest[chain]
                    for subchain in subchains:
                        single_chains = list(subchain)
                        for single in single_chains:
                            outdict = self.retrieve_residues(single, range_of_interest[cdr])
                            pdb_chunks[(chain, cdr)] = outdict  # Store the results in a dictionary
        print(f"PDB lines retrieved for {self.name}: {len(pdb_chunks)}")
        return pdb_chunks  # Return the dictionary of PDB lines

    def get_coords(self, line):
        x = line[30:38].strip()
        y = line[38:46].strip()
        z = line[46:54].strip()
        return np.array([float(x), float(y), float(z)])

    def dist(self, residue_one, residue_two):
        diff = residue_one - residue_two  # Calculate the difference between two vectors
        return np.sqrt(np.sum(diff * diff))  # Return the Euclidean distance

    def get_atom(self, line):
        return line[12:16].strip()  # Return the atom name

    def get_residue(self, line):
        return line[22:26].strip()  # Return the residue sequence number

    def get_aa(self, line):
        three2one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                     'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        AA = line[17:20].strip()  # Extract three-letter amino acid
        return three2one[AA]  # Convert three-letter code to one-letter code

    def calculate_distance(self, from_list, to_list):
        data = []

        for row, i in enumerate(from_list):
            from_vector = self.get_coords(i)  # Get coordinates for the first list
            from_atom = self.get_atom(i)  # Get atom name
            from_residue = self.get_residue(i)  # Get residue number
            from_aa = self.get_aa(i)  # Get amino acid

            for col, j in enumerate(to_list):
                to_vector = self.get_coords(j)  # Get coordinates from the first list
                to_atom = self.get_atom(j)  # Get atom name
                to_residue = self.get_residue(j)  # Get residue number
                to_aa = self.get_aa(j)  # Get amino acid

                distance = self.dist(from_vector, to_vector)  # Calculate distance between two points

                if distance <= 5:  # If distance is within 5 angstroms
                    data.append([from_atom, from_residue, from_aa,
                                 to_atom, to_residue, to_aa,
                                 distance])  # Append data to the list
        print(f"Distances calculated: {len(data)}")
        return pd.DataFrame(data, columns=["donor_atom", "donor_residue", "donor_aa",
                                           "to_atom", "to_residue", "to_aa",
                                           "distance"])

    def calculate_all_distances(self):
        if self.too_many_atoms:
            print("PDB file is too large, skipping")
            return None

        matrices = []

        for i in self.immuno_pdb:
            if not self.immuno_pdb[i]:  # Check if list is empty
                print(f"No immunoglobulin residues found for chain {i} in {self.name}")
                continue
            for j in self.antigen_pdb:
                if not self.antigen_pdb[j]:  # Check if list is empty
                    print(f"No antigen residues found for chain {j} in {self.name}")
                    continue
                matrix = self.calculate_distance(self.immuno_pdb[i], self.antigen_pdb[j])
                if not matrix.empty:
                    matrix["donor_chain"] = i[0]
                    matrix["donor_loop"] = i[1]
                    matrix["acceptor_chain"] = j[0]
                    matrix["acceptor_loop"] = j[1]
                    matrix["Pdb"] = self.name
                    matrices.append(matrix)

        if matrices:
            self.distances = pd.concat(matrices)
            print(f"Distances matrix created for {self.name}")
        else:
            print(f"No distances found for {self.name}")

    def write(self):
        if self.distances is not None:
            print(f"Writing to {self.outfile}")
            self.distances.to_csv(self.outfile, index=False)
        else:
            print(f"No distances calculated for {self.path}")

def caller(immuno_chains, antigen_chains, pdb_file):
    outfile = pdb_file.replace(".pdb", ".csv")
    outfile = outfile.replace("imgt", "csvs")
    ab = Immunoglobulin(pdb_file, immuno_chains, antigen_chains, outfile)
    ab.calculate_all_distances()
    ab.write()

def main():
    n_cores = mp.cpu_count()
    print("Number of processors:", n_cores)

    # Define paths to PDB files
    pdb_directory = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/imgt/*.pdb"
    pdb_files = glob.glob(pdb_directory)

    # Ensure output directory exists
    output_directory = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/csvs"
    os.makedirs(output_directory, exist_ok=True)

    antibody_chains = {"heavy": "HCHAIN",
                       "light": "LCHAIN"}

    antigen_chains = {"antigen": "AGCHAIN"}

    print(f"Analyzing {len(pdb_files)} PDB files.")

    Parallel(n_jobs=n_cores, verbose=5)(
        delayed(caller)(antibody_chains, antigen_chains, pdb_file)
        for pdb_file in pdb_files
    )

if __name__ == "__main__":
    main()
