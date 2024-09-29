import simpy
import random
import numpy as np
import matplotlib.pyplot as plt

# Simulate a single RNA strand evolution with simple mutation and fitness evaluation
def rna_evolution(env, generations, mutation_rate, target_sequence):
    rna_sequence = ''.join(random.choices('ACGU', k=len(target_sequence)))
    for _ in range(generations):
        # Mutate the sequence
        rna_sequence = ''.join(
            random.choice('ACGU'.replace(base, '')) if random.random() < mutation_rate else base
            for base in rna_sequence
        )
        yield env.timeout(1)

    # Calculate fitness at the end of the evolution
    return sum(a == b for a, b in zip(rna_sequence, target_sequence)) / len(target_sequence)


# Run a single simulation and return the fitness score
def run_simulation(generations, mutation_rate, target_sequence):
    env = simpy.Environment()
    fitness = env.process(rna_evolution(env, generations, mutation_rate, target_sequence))
    env.run()
    return fitness.value


if __name__ == "__main__":
    # Parameters for the hackathon
    target_sequence = "ACGUACGUACGU"
    generations = 100
    mutation_rate = 10**-6

    # Monte Carlo Simulation
    N = 10000  # Reduced number of simulations for a quick demo
    fitness_results = [run_simulation(generations, mutation_rate, target_sequence) for _ in range(N)]

    # Quick analysis and display
    mean_fitness = np.mean(fitness_results)
    std_fitness = np.std(fitness_results)

    print(f'Mean Fitness: {mean_fitness}, Standard Deviation: {std_fitness}')

    plt.hist(fitness_results, bins=10)
    plt.title('Fitness Distribution after RNA Evolution')
    plt.xlabel('Fitness')
    plt.ylabel('Frequency')
    plt.show()
