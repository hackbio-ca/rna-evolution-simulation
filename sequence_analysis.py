########## Sequence Length ###########

import matplotlib.pyplot as plt

generations = 5 # 
sequences = ['''tgttttaca ggttcgcgac gtgctcgtac gtggctttgg
      agactccgtg gaggaggtct tatcagaggc acgtcaacat cttaaagatg gcacttgtgg
      cttagtagaa gttgaaaaag gcgttttgcc tcaacttgaa cagccctatg tgttcatcaa
      acgttcggat gctcgaactg cacctcatgg tcatgttatg gttgagctgg tagcagaact
      cgaaggcatt cagtacggtc gtagtggtga gacacttggt gtccttgtcc ctcatgtggg
       cgaaatacca gtggcttacc gcaaggttct tcttcgtaag aacggtaata aaggagctgg
      tggccatagt tacggcgccg atctaaagtc atttgactta ggcgacgagc ttggcactga
      tccttatgaa gattttcaag aaaactggaa cactaaacat agcagtggtg''',
      '''tgttttaca ggttcgcgac gtgctcgtac gtggctttgg
      agactccgtg gaggaggtct tatcagaggc acgtcaacat cttaaagatg gcacttgtgg
      cttagtagaa gttgaaaaag ccttgtcc ctcatgtggg
       cgaaatacca gtggcttacc gcaaggttct tcttcgtaag aacggtaata aaggagctgg
      tggccatagt tacggcgccg atctaaagtc atttgactta ggcgacgagc ttggcactga
      tccttatgaa gattttcaag aaaactggaa cactaaacat agcagtggtg''',
      '''tgttttaca ggttcgcgac gtgctcgtac gtggctttgg
      agactccgtg gaggaggtct tatcagaggc acgtcaacat cttaaagatg gcacttgtgg
      cttagtagaa gttgaaaaag ccttgtcc ctcatgtggg
       cgaaatacca gtggcttacc gcaaggttct tcttcgtaag aacggtaata aaggagctgg
      tggccatagt tacggcgccg atctaaagtc atttgactta ggcgacgagc ttggcactga
      tccttatgaa gattttcaag aaaactggaa cactaaacat agcagtggtg''',
      '''tgttttaca ggttcgcgac gtgctcgtac gtggctttgg
      agactccgtg gaggaggtct tatcagaggc acgtcaacat cttaaagatg gcacttgtgg
      cttagtagaa gttgaaaaag ccttgtcc ctcatgtggg
       cgaaatacca gtggcttacc gcaaggttct tcttcgtaag aacggtaata aaggagctgg
      tggccatagt tacggcgccg atctaaagtc atttgactta ggcgacgagc ttggcactga
      tccttatgaa gattttcaag aaaactggaa cactaaacat agcagtggtg''',
      '''tgttttaca ggttcgcgac gtgctcgtac gtggctttgg
      agactccgtg gaggaggtct tatcagaggc acgtcaacat cttaaagatg gcacttgtgg
      cttagtagaa gttgaaaaag ccttgtcc ctcatgtggg
     '''
      ]


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
query = 'UUCGG'

aligner = Align.PairwiseAligner()
# aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

alignments = aligner.align(original, query)

if len(original) >= len(query):
    max = aligner.align(original, original).score
else:
    max = aligner.align(query, query).score

similarity_percentage = alignments.score / max * 100
print(similarity_percentage)

# Visualizing the alignment
# for alignment in sorted(alignments):
#     print("Score = %.1f:" % alignment.score)
#     print(alignment)


######### Function Analysis ##########

from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
import threading
import time


# for i in range(generations): 
#     # sequence = '> Sequence ' + str(i + 1) + sequences[i]
#     sequence = Seq(sequences[i])

#     result = NCBIWWW.qblast("blastx", "nr", sequence) 

#     E_VALUE_THRESH = 1e-20 
#     for record in NCBIXML.parse(result): 
#         if record.alignments: 
#             print("\n") 
#             print("query: %s" % record.query[:100]) 
#             found = False
#             for align in record.alignments: 
#                 for hsp in align.hsps: 
#                     if hsp.expect < E_VALUE_THRESH: 
#                         print("match: %s " % align.title[:100])
#                         # aligned_sequence = align.title[:100]
#                         found = True
#                         break
#                 if found:
#                     break

threads = [None] * generations
results = [None] * generations

def process_blast(results, index, *args):
    print('processing sequence: ' + str(args[2]) + 'at ' + 'index ' + str(index))
    result = NCBIWWW.qblast(*args)
    results[index] = result

start = time.time()

for i in range(generations):
    # sequence = '> Sequence ' + str(i + 1) + sequences[i]
    sequence = Seq(sequences[i])
    threads[i] = threading.Thread(target=process_blast, args=(results, i, "blastx", "nr", sequence))
    threads[i].start()

for thread in threads:
    thread.join()

print(time.time() - start)

for result in results:
    E_VALUE_THRESH = 1e-20
    for record in NCBIXML.parse(result):
        if record.alignments:
            print("\n") 
            print("query: %s" % record.query[:100])
            found = False
            for align in record.alignments:
                for hsp in align.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        print("match: %s " % align.title[:100])
                        # aligned_sequence = align.title[:100]
                        found = True
                        # Comment all the code below if you want to get 10 BLAST results
                        break
                if found:
                    break

