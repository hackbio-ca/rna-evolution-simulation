import math

import simpy
import random

# RNA nucleotides
nucleotides = ['A', 'U', 'C', 'G']


# Function to simulate RNA replication with mutation
def replicate_rna(env, initial_rna, mutation_rate, generations):
    rna_population = [initial_rna]

    # Loop over each generation
    for gen in range(generations):
        print(f"--- Generation {gen + 1} ---")
        new_population = []

        # Replicate each RNA in the current population
        for rna in rna_population:
            new_rna = ""
            for nucleotide in rna:
                # Apply mutation with a certain probability
                if random.random() < mutation_rate:
                    mutated_nucleotide = random.choice(nucleotides)
                    new_rna += mutated_nucleotide
                    print(f"Mutation: {nucleotide} -> {mutated_nucleotide}")
                else:
                    new_rna += nucleotide

            print(f"Original RNA: {rna}")
            print(f"New RNA: {new_rna}")
            new_population.append(new_rna)

        # Update the population for the next generation
        rna_population = new_population

        # Simulate the time until the next generation (optional, can be modified)
        yield env.timeout(1)


if __name__ == "__main__":
    # Simulation environment
    env = simpy.Environment()

    # Parameters for the replication
    initial_rna = "AUGCGUAAU"  # Starting RNA sequence
    mutation_rate = 0.01      # Mutation rate (10% per nucleotide)
    generations = 20            # Number of replication generations

    # Start the RNA replication process
    env.process(replicate_rna(env, initial_rna, mutation_rate, generations))

    # Run the simulation
    env.run()
