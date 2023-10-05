"""Test cli.py functions."""
from pathlib import Path
import pytest  # type: ignore
from scRNAsim_toolz.fragment_selector.cli import (
    file_validation, main
)

TEST_FILES_DIR = Path(__file__).resolve().parent / "test_files"
FASTA_FILE = TEST_FILES_DIR / "test.fasta"
CSV_FILE = TEST_FILES_DIR / "test.csv"
TAB_FILE = TEST_FILES_DIR / "test_tab.csv"

_, csv_counts = file_validation(str(FASTA_FILE), str(CSV_FILE), ",")
_, tab_counts = file_validation(str(FASTA_FILE), str(TAB_FILE), "\t")


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
    with pytest.raises(SystemExit):
        main()
