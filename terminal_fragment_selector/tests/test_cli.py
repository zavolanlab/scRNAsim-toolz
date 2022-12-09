"""Test cli.py functions."""
import pytest
# from Bio import SeqIO

from term_frag_sel.cli import file_validation  # type: ignore

FASTA_FILE = "tests/test_files/test.fasta"


def test_file():
    """Test check_positive function."""
    with pytest.raises(FileNotFoundError):
        file_validation("", "", ",")
