import glob as g
import multiprocessing as mp
import os
from joblib import Parallel, delayed

class Immunoglobulin:
    total = 0

    def __init__(self, path, im_chains, ag_chains):
        Immunoglobulin.total += 1
        self.path = path
        self.name = os.path.basename(path).replace(".pdb", "")
        self.sc_file = self.name + "_scpipe.txt"
        self.sc_out = self.name + "_scout.txt"
        self.immuno_chains = self.get_chains(im_chains)
        self.antigen_chains = self.get_chains(ag_chains)

    def get_chains(self, chains):
        data = {}
        try:
            with open(self.path) as f:
                for line in f:
                    if line.startswith("REMARK") and "=" in line:
                        line = line.split()
                        for chain in chains:
                            chainstr = chains[chain]
                            index = list(filter(lambda x: x.startswith(chainstr), line))
                            if len(index) != 0:
                                residues = index[0]
                                residues = residues.split("=")[-1]
                                if residues != "NONE":
                                    residues = residues.split(";")
                                    data[chain] = residues
        except Exception as e:
            print(f"Error reading chains from {self.path}: {e}")
        return data

    def write_sc_pipe(self):
        try:
            with open(self.sc_file, "w") as infile:
                infile.write("molecule 1\n")
                for chain in self.immuno_chains:
                    chain = self.immuno_chains[chain]
                    for subchain in chain:
                        for individual in list(subchain):
                            infile.write("chain %s\n" % individual)

                infile.write("molecule 2\n")
                for chain in self.antigen_chains:
                    chain = self.antigen_chains[chain]
                    for subchain in chain:
                        for individual in list(subchain):
                            infile.write("chain %s\n" % individual)

                infile.write("end\n")
        except Exception as e:
            print(f"Error writing SC pipe file {self.sc_file}: {e}")

    def run_SC(self):
        SC = None
        try:
            cmdstring = "sc XYZIN %s < %s > %s" % (self.path, self.sc_file, self.sc_out)
            os.system(cmdstring)

            with open(self.sc_out, 'r') as f:
                for line in f:
                    if 'Shape complementarity statistic Sc' in line:
                        SC = line.strip().split("=")[-1]
        except Exception as e:
            print(f"Error running SC for {self.path}: {e}")
        return SC

    def cleanup(self):
        try:
            os.remove(self.sc_file)
            os.remove(self.sc_out)
        except Exception as e:
            print(f"Error cleaning up files {self.sc_file} and {self.sc_out}: {e}")

def process_antibody(antibody, antibody_chains, antigen_chains, outfile_path):
    try:
        print(f"Processing antibody: {antibody}")
        ab = Immunoglobulin(antibody, antibody_chains, antigen_chains)

        if ab.immuno_chains and ab.antigen_chains:
            ab.write_sc_pipe()
            sc = ab.run_SC()
            ab.cleanup()
            with open(outfile_path, "a") as outfile:
                outfile.write(ab.name + "," + str(sc) + ",antibody\n")
            print(f"Completed: {ab.name}, SC: {sc}")
        else:
            print(f"Chains not found for {ab.name}")
    except Exception as e:
        print(f"Error processing {antibody}: {e}")

def main():
    n_cores = mp.cpu_count()
    print("Number of processors:", n_cores)

    antibodies = g.glob("/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/imgt/*.pdb")
    no_antibodies = len(antibodies)
    print(f"There are {no_antibodies} antibodies")

    antibody_chains = {"heavy": "HCHAIN", "light": "LCHAIN"}
    antigen_chains = {"antigen": "AGCHAIN"}

    outfile_path = "/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/sc_scores.csv"
    try:
        with open(outfile_path, "w") as outfile:
            outfile.write("name,sc,origin\n")
    except Exception as e:
        print(f"Error opening output file {outfile_path}: {e}")

    Parallel(n_jobs=n_cores, backend='threading', verbose=5)(
        delayed(process_antibody)(antibody, antibody_chains, antigen_chains, outfile_path) for antibody in antibodies
    )

if __name__ == "__main__":
    main()
