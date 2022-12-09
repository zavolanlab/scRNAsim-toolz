"""Test utils.py functions."""
import pytest
from Bio import SeqIO

from term_frag_sel.fragmentation import fragmentation  # type: ignore

with open("tests/test_files/test.fasta", "r", encoding="utf-8") as handle:
    fasta = SeqIO.parse(handle, "fasta")

# have to create the counts file

MU = 100
STD = 10
A_PROB = 0.22
C_PROB = 0.28
T_PROB = 0.25
G_PROB = 0.25


def test_frag():
    """Test fragmentation function."""
    with pytest.raises(TypeError):
        fragmentation()
