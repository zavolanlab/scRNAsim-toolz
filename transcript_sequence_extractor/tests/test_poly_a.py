import pytest
from poly_a import poly_a_generator, poly_a_addition_to_fasta_list




def test_poly_a_generator():
    assert poly_a_generator(exon_string) == exon_and_polya


def test_poly_a_addition_to_fasta_list():
    assert poly_a_addition_to_fasta_list(list_of_tuples) == manipulated_list_of_tuples
