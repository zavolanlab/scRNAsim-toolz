import numpy as np

def PolyA_generator(exon):
def PolyA_generator(exon):
"""Appends PolyA tail.

    This function appends a polyA tail on the transcript sequence based on a random selection of bases, with A defined as being in larger proportion..

    Parameters
    ----------
    arg1 : exon
             Transcript sequence made up of exons.

    Returns
    -------
    fasta
        Transcript sequence with concatenated polyA tails.



    Raises
    ------
    TypeError
        ValueError: PolyA tail needs to contain majority As.
    """
listA = ['A','T','G','C']
        polyA = ''.join(np.random.choice(listA,250,p=[0.9,0.040,0.020,0.020]))
        return (exon+polyA)
	listA = ['A','T','G','C']
	polyA = ''.join(np.random.choice(listA,250,p=[0.9,0.040,0.020,0.020]))
	return (exon+polyA)
