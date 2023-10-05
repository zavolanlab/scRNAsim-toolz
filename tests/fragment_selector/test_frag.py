"""Test utils.py functions."""
from pathlib import Path
import pandas as pd  # type: ignore
import pytest  # type: ignore
from Bio import SeqIO  # type: ignore
from scRNAsim_toolz.fragment_selector.fragmentation import (
    fragmentation
)

TEST_FILES_DIR = Path(__file__).resolve().parent / "test_files"
FASTA_FILE = TEST_FILES_DIR / "test.fasta"
CSV_FILE = TEST_FILES_DIR / "test.csv"

fasta_dict = {}
with open(FASTA_FILE, "r", encoding="utf-8") as handle:
    fasta_sequences = SeqIO.parse(handle, "fasta")
    for record in fasta_sequences:
        fasta_dict[record.id] = str(record.seq).upper()

seq_counts = pd.read_csv(CSV_FILE, names=["seqID", "count"])
seq_counts = seq_counts.astype({"seqID": str})


def test_fragmentation():
    """Test fragmentation function."""
    assert isinstance(fragmentation(fasta_dict, seq_counts), list)

    with pytest.raises(ValueError):
        fragmentation(fasta_dict, seq_counts, mean='a')

    with pytest.raises(ValueError):
        fragmentation(fasta_dict, seq_counts, std='a')
