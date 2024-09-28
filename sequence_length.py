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

