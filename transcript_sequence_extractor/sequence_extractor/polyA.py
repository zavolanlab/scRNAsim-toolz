import numpy as np
# To do: Taking probabilities of nucleotides from user and raising error if sum != 1
def PolyA_generator(
	exon: str,
) -> str:
	"""Adds a PolyA tail to an exon sequence input into the function.

	 Args:
		exon: RNA sequence, obtained from concatenation of exons, that needs polyA to be added to its 3' end.

	Returns:
		RNA with polyA tail added to its 3' end.
	"""
	listA = ['A','T','G','C']
	polyA = ''.join(np.random.choice(listA,250,p=[0.9,0.040,0.020,0.020]))
	return (exon+polyA)
