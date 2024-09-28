########## Sequence Length ###########

import matplotlib.pyplot as plt

generations = 5
sequences = ['ATTGC', 'ATTGC', 'ATTCG', 'ATCG', 'ATCCG']

x = []
y = []

for i in range(1, generations + 1):
    x.append(i)

for item in sequences:
    y.append(len(item))

plt.plot(x, y) 
plt.xlabel("Generations")  
plt.ylabel("Sequence Length") 
plt.title("Sequence Length Throughout Generations of Mutations") 
plt.show()

######### Sequence Similarity ##########
from Bio import Align
from Bio.Align import substitution_matrices

original = 'UUGG'
query = 'GCU'

aligner = Align.PairwiseAligner()
# aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

alignments = aligner.align(original, query)

if len(original) >= len(query):
    max = aligner.align(original, original).score
else:
    max = aligner.align(query, query).score

similarity_percentage = alignments.score / max * 100
print(similarity_percentage)

# for alignment in sorted(alignments):
#     print("Score = %.1f:" % alignment.score)
#     print(alignment)

