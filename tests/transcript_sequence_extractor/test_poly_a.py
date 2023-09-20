"""Test poly_a.py script."""

import numpy as np
from sequence_extractor.poly_a import (
    poly_a_generator, poly_a_addition_to_fasta_list
)


def test_poly_a_generator():
    """Test poly_a_generator function."""
    exon_string = 'ACGGCTGGCACCTTGTTTGGGGAAGGATTCCGTGCCTTTG'
    poly_a_length = np.random.randint(100, 250)
    test_string = poly_a_generator(exon_string, poly_a_length)

    assert len(test_string) == sum([len(exon_string), poly_a_length])


def test_poly_a_addition_to_fasta_list():
    """Test poly_a_addition_to_fasta_list function."""
    exon_list = [('>ENST00000673477', 'ACGGCTGGCACCTTGTTTGGGGAAG'),
                 ('>ENST00000378391', 'AAATACTGACGGACGTGG')]
    poly_a_length = np.random.randint(100, 250)
    test_list = poly_a_addition_to_fasta_list(exon_list, poly_a_length)

    assert len(test_list[0][1]) == sum([len(exon_list[0][1]), poly_a_length])
    assert len(test_list[1][1]) == sum([len(exon_list[1][1]), poly_a_length])
