import numpy as np

def PolyA_generator(exon):
	listA = ['A','T','G','C']
	polyA = ''.join(np.random.choice(listA,250,p=[0.9,0.040,0.020,0.020]))
	return (exon+polyA)
