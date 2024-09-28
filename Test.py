import math

import simpy
import random

# RNA nucleotides
nucleotides = ['A', 'U', 'C', 'G']
unique = []

# Function to simulate RNA replication with mutation
def replicate_rna(env, initial_rna, mutation_rate, generations):
    rna_population = [initial_rna]

    # Loop over each generation
    for gen in range(generations):
        #print(f"--- Generation {gen + 1} ---")
        new_population = []

        # Replicate each RNA in the current population
        for rna in rna_population:
            new_rna = ""
            for nucleotide in rna:
                # Apply mutation with a certain probability
                if random.random() < mutation_rate:
                    mutated_nucleotide = random.choice(nucleotides)
                    new_rna += mutated_nucleotide
                    #print(f"Mutation: {nucleotide} -> {mutated_nucleotide}")
                else:
                    new_rna += nucleotide

            #print(f"Original RNA: {rna}")
            #print(f"New RNA: {new_rna}")
            new_population.append(new_rna)

            if new_rna not in unique:
                unique.append(new_rna)

        # Update the population for the next generation
        rna_population = new_population

        # Simulate the time until the next generation (optional, can be modified)
        yield env.timeout(1)

    entropy()

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import subprocess

def entropy():
    sequences = []
    for seq in unique:
        seq_rec = SeqRecord(Seq(seq), id=seq)
        sequences.append(seq_rec)

    with open("rna_sequences.fasta", "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

    # Specify the input FASTA file and output alignment file
    in_file = "rna_sequences.fasta"
    out_file = "aligned_rna_sequences.aln"

    # Create a command line object for Clustal Omega
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)

    # Run Clustal Omega
    stdout, stderr = clustalomega_cline()

    # Read the alignment from the output file
    alignment = AlignIO.read(out_file, "clustal")

    # Print the alignment
    print(alignment)




if __name__ == "__main__":
    # Simulation environment
    env = simpy.Environment()

    # Parameters for the replication
    initial_rna = "AUGCGUAAU"  # Starting RNA sequence
    mutation_rate = 10**-4      # Mutation rate (10% per nucleotide)
    generations = 1000            # Number of replication generations

    # Start the RNA replication process
    env.process(replicate_rna(env, initial_rna, mutation_rate, generations))

    # Run the simulation
    env.run()
