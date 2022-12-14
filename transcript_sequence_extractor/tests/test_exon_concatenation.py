import pytest
import exon_concatenation from exon_concatenation

test_fasta_1 = "test_files/test_1.fa"
test_fasta_2 = "test_files/test_2.fa"

def test_exon_concatenation():
    assert exon_concatenation(test_fasta_1) == expected_list_of_tuples
    assert exon_concatenation(test_fasta_2) == expected
