import os
import pandas as pd
import glob as g
from joblib import Parallel, delayed

class Immunoglobulin:
    total = 0

    def __init__(self, path, im_chains, ag_chains, outfile):
        Immunoglobulin.total += 1
        self.path = path
        self.name = os.path.basename(path).replace(".pdb", "")
        self.outfile = outfile
        self.immuno_chains = self.get_chains(im_chains)
        self.antigen_chains = self.get_chains(ag_chains)
        self.cdr_ranges = {"CDR1": range(27, 39), "CDR2": range(56, 66), "CDR3": range(105, 118)}
        self.antigen_ranges = {"antigen": range(0, 10000)}
        self.too_many_atoms = self.count_atoms()
        self.immuno_pdb = self.retrieve_pdb_lines(self.immuno_chains, self.cdr_ranges)
        self.antigen_pdb = self.retrieve_pdb_lines(self.antigen_chains, self.antigen_ranges)
        self.kmers = []

        print(f"Initialized Immunoglobulin object for {self.name}")
        print(f"Immuno chains: {self.immuno_chains}")
        print(f"Antigen chains: {self.antigen_chains}")

    def count_atoms(self):
        atom_count = 0
        with open(self.path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    atom_count += 1
                    if atom_count > 100000:
                        print(f"PDB file {self.name} has more than 100000 atoms.")
                        return True
        print(f"PDB file {self.name} has {atom_count} atoms.")
        return False

    def get_chains(self, chains):
        data = {}
        with open(self.path) as f:
            for line in f:
                if line.startswith("REMARK") and "=" in line:
                    line = line.split()
                    for chain in chains:
                        chainstr = chains[chain]
                        index = list(filter(lambda x: x.startswith(chainstr), line))
                        if index:
                            residues = index[0].split("=")[-1]
                            if residues != "NONE":
                                residues = residues.split(";")
                                data[chain] = residues
        print(f"Retrieved chains: {data}")
        return data

    def retrieve_residues(self, c, res):
        lines = []
        amino_acids = ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 'ILE', 'PRO', 'THR', 'PHE', 'ASN', 'ALA', 'VAL', 'GLU', 'TYR', 'MET']
        with open(self.path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain = line[21].strip() == c
                    residues = int(line[22:26].strip()) in res
                    amino_acid = line[17:20].strip() in amino_acids
                    if chain and residues and amino_acid:
                        lines.append(line)
        print(f"Retrieved residues for chain {c} in {res}: {lines[:5]}...")
        return lines

    def retrieve_pdb_lines(self, chains_of_interest, range_of_interest):
        if self.too_many_atoms:
            print(f"PDB file {self.name} is too large, skipping")
            return None
        pdb_chunks = {}
        if chains_of_interest:
            for chain in chains_of_interest:
                for cdr in range_of_interest:
                    subchains = chains_of_interest[chain]
                    for subchain in subchains:
                        single_chains = list(subchain)
                        for single in single_chains:
                            outdict = self.retrieve_residues(single, range_of_interest[cdr])
                            pdb_chunks[(chain, cdr)] = outdict
        print(f"Retrieved PDB lines: {pdb_chunks}")
        return pdb_chunks

    def create_3mers(self):
        amino_acids = ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 'ILE', 'PRO', 'THR', 'PHE', 'ASN', 'ALA', 'VAL', 'GLU', 'TYR', 'MET']
        three2one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        if self.antigen_pdb:
            for chain, residues in self.antigen_pdb.items():
                sequence = []
                for line in residues:
                    if line.startswith("ATOM"):
                        amino_acid = line[17:20].strip()
                        if amino_acid in amino_acids:
                            sequence.append(three2one[amino_acid])
                sequence = "".join(sequence)
                for pos in range(len(sequence) - 2):  # Adjust window for 3-mers
                    kmer = sequence[pos:pos + 3]
                    self.kmers.append((chain[0], pos, kmer))

    def create_data_frame(self):
        headers = ["chain_id", "position", "kmer"]
        kmer_df = pd.DataFrame(self.kmers, columns=headers)
        os.makedirs(os.path.dirname(self.outfile), exist_ok=True)
        print(f"Number of k-mers: {len(self.kmers)}")
        if not self.kmers:
            print("No k-mers found. Exiting.")
            return
        print(f"Writing k-mers to {self.outfile}")
        kmer_df.to_csv(self.outfile, index=False)
        print(f"Successfully wrote k-mers to {self.outfile}")

def caller(immuno_chains, antigen_chains, pdb_file):
    outfile = pdb_file.replace(".pdb", "_3mers.csv").replace("imgt", "kmers")
    try:
        ab = Immunoglobulin(pdb_file, immuno_chains, antigen_chains, outfile)
        ab.create_3mers()
        ab.create_data_frame()
    except Exception as e:
        print(f"Error processing file {pdb_file}: {e}")

def main():
    n_cores = 64
    print(f"Number of processors: {n_cores}")
    pdb_directory = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/imgt/*.pdb"
    antibodies = g.glob(pdb_directory)
    no_antibodies = len(antibodies)
    print(f"There are {no_antibodies} antibodies")
    antibody_chains = {"heavy": "HCHAIN", "light": "LCHAIN"}
    antigen_chains = {"antigen": "AGCHAIN"}
    print("Analysing antibodies")
    Parallel(n_jobs=n_cores, backend='threading', verbose=5)(delayed(caller)(antibody_chains, antigen_chains, pdb_file) for pdb_file in antibodies)

if __name__ == "__main__":
    main()
