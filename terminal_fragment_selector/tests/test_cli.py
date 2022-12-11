"""Test cli.py functions."""
import pytest

from term_frag_sel.cli import file_validation, main  # type: ignore

FASTA_FILE = "tests/test_files/test.fasta"
CSV_FILE = "tests/test_files/test.csv"
TAB_FILE = "tests/test_files/test_tab.csv"

_, csv_counts = file_validation(FASTA_FILE, CSV_FILE, ",")
_, tab_counts = file_validation(FASTA_FILE, TAB_FILE, "\t")


def test_file():
    """Test file_validation function."""
    assert isinstance(file_validation(FASTA_FILE, CSV_FILE, ","), tuple)
    assert isinstance(file_validation(FASTA_FILE, TAB_FILE, "\t"), tuple)
    assert csv_counts.equals(tab_counts)

    with pytest.raises(FileNotFoundError):
        file_validation(FASTA_FILE, "", ",")
    with pytest.raises(FileNotFoundError):
        file_validation("", CSV_FILE, ",")
    with pytest.raises(ValueError):
        file_validation(CSV_FILE, CSV_FILE, ",")


def test_main():
    """Test main() function."""
    with pytest.raises(TypeError):
        main("")
