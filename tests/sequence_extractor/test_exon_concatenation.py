"""Test exon_concatenation.py."""
from pathlib import Path
from ...scRNAsim_toolz.sequence_extractor.exon_concat import (  # type: ignore
    exon_concatenation
)

test_dir = Path(__file__).parent.resolve()

test_fasta_1 = test_dir / "test_files" / "test_1.fa"
test_fasta_2 = test_dir / "test_files" / "test_2.fa"


def test_exon_concatenation():
    """Test exon_concatenation function."""
    # Test for test_fasta_1
    expected_fasta_1 = [
        (">ENST00000673477",
         "TTTCGCCTGCGCAGTGGTCCTGGCCACCGGCTCGCGGCGCGTGGAGGCTGC"
         "TCCCAGCCGCGCCCGAGTCAGACTCGGGTGGGGGTCCCGGTTACGCCAAGG"
         "AGGCCCTGAATCTGGCGCAGATGCAGGAGCAGACGCTGCAGTTGGAGCAAC"
         "AGTCCAAGCTCAAA"),
        (">ENST00000378391",
         "AAATACTGACGGACGTGGAAGTGTCGCCCCAGGAAGGCTGCATCACAAAGTC"
         "TCCGAAGACCTGGGCAGTGAGAAGTTCTGCGTGGATGCAAATCAGGCGGGGG"),
    ]

    result_fasta_1 = exon_concatenation(test_fasta_1)
    assert result_fasta_1 == expected_fasta_1

    # Test for test_fasta_2
    expected_fasta_2 = [
        (">ENST00000673477",
         "ACGGCTGGCACCTTGTTTGGGGAAGGATTCCGTGCCTTTGTGACAGACCGGGA"
         "CAAAGTGACTGGCTGGGCTGACGCTGCTGGCTGTCGGGGTCTACTCAGCCAAG"
         "AATGCGATCAGCCGGCGGCTCCTCAGTCGACCCCAGGACGTGCTGGAGGGTGT"
         "TGTGCTTAGT"),
    ]

    result_fasta_2 = exon_concatenation(test_fasta_2)
    assert result_fasta_2 == expected_fasta_2
