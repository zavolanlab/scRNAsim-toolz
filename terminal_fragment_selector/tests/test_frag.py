"""Test utils.py functions."""
import pandas as pd
import pytest
from Bio import SeqIO

from term_frag_sel.fragmentation import fragmentation

FASTA_FILE = "tests/test_files/test.fasta"
CSV_FILE = "tests/test_files/test.csv"

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
