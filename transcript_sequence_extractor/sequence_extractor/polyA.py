import numpy as np
# To do: Taking probabilities of nucleotides from user and raising error if sum != 1
def polyA_generator(
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

def polyA_addition_to_fasta_list(
	fasta_list: list,
) -> list:
	"""Takes in a list of alternate annotations and exons and outputs a list where polyA tail has been added to all the exon 3' ends.

	Args:
		fasta_list: List contaning annotations and exons, with every element being an equivalent to a line in the fasta file

	Returns:
		A list like the initial list, this time with polyA tail added onto it.
	"""
	mature_rna_list = [(i[0],polyA_generator(i[1])) for i in fasta_list]
	return mature_rna_list
