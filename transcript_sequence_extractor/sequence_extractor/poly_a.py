""" This script contains two functions and the first function is called by the second function and used to add poly A tail to the concatenated exon"""
import numpy as np
# To do: Taking probabilities of nucleotides from user and raising error if sum != 1
def poly_a_generator(
	exon: str,
) -> str:
    """Adds a PolyA tail to an exon sequence input into the function.

	 Args:
		exon: RNA sequence, obtained from concatenation of exons, that needs polyA to be added to its 3' end.

	Returns:
		RNA with polyA tail added to its 3' end.
	"""
    list_of_nucleotides = ['A','T','G','C']
    poly_a_string = ''.join(np.random.choice(list_of_nucleotides,250,p=[0.914,0.028,0.025,0.033]))
    return exon+poly_a_string

def poly_a_addition_to_fasta_list(
	fasta_list: list,
) -> list:
    """Takes in a list of tuples with annotations and exons and outputs a list where polyA tail has been added to all the exon 3' ends.

	Args:
		fasta_list: List contaning tuples of annotations and exons

	Returns:
		A list like the initial list, this time with polyA tail added onto it.
	"""
    mature_rna_list = [(i[0],poly_a_generator(i[1])) for i in fasta_list]
    return mature_rna_list
