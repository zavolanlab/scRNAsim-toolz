fa = open("fasta.fa",'r')
lines = fa.readlines()
for x in range(int(len(lines)/2)):
    if x == 0:
        annotation = lines[0]
        read = lines[1]
    if x >= 1:
        if lines[2*x] == lines[2*(x-1)]:
            read+= lines[(2*x)+1]
        else:
            annotation = lines[2*x]
            read = lines[(2*x)+1]

# Function for random addition of polyA to sequences

import numpy as np

listA = ['A','U','G','C']
''.join(np.random.choice(listA,250,p=[0.9,0.040,0.020,0.020]))
