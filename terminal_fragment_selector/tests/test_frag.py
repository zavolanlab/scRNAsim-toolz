"""Test utils.py functions."""
import pandas as pd
import pytest
from Bio import SeqIO

from term_frag_sel.fragmentation import fragmentation, get_cut_number

FASTA_FILE = "tests/test_files/test.fasta"
CSV_FILE = "tests/test_files/test.csv"

fasta_dict = {}
with open(FASTA_FILE, "r", encoding="utf-8") as handle:
    fasta_sequences = SeqIO.parse(handle, "fasta")
    for record in fasta_sequences:
        fasta_dict[record.id] = str(record.seq).upper()

seq_counts = pd.read_csv(CSV_FILE, names=["seqID", "count"])
seq_counts = seq_counts.astype({"seqID": str})

NUC = {'A': 0.22, 'T': 0.25, 'G': 0.25, 'C': 0.28}


def test_get_cut_number():
    """Test get_cut_number function."""
    assert get_cut_number(100, 20) <= 15 and get_cut_number(100, 20) >= 0

    with pytest.raises(ValueError):
        get_cut_number(seq_len='a', mean=4)
    with pytest.raises(ValueError):
        get_cut_number(seq_len=10, mean='a')


def test_fragmentation():
    """Test fragmentation function."""
    assert isinstance(fragmentation(fasta_dict, seq_counts,
                                    NUC), list)

    # no need to check string mean or std since it's checked at CLI
    with pytest.raises(ValueError):
        nuc_probs = {'A': 'a', 'T': 0.25,
                     'G': 0.25, 'C': 0.28}
        fragmentation(fasta_dict, seq_counts, nuc_probs)
    with pytest.raises(ValueError):
        nuc_probs = {'A': 0.22, 'T': 'a',
                     'G': 0.25, 'C': 0.28}
        fragmentation(fasta_dict, seq_counts, nuc_probs)
    with pytest.raises(ValueError):
        nuc_probs = {'A': 0.22, 'T': 0.25,
                     'G': 'a', 'C': 0.28}
        fragmentation(fasta_dict, seq_counts, nuc_probs)
    with pytest.raises(ValueError):
        nuc_probs = {'A': 0.22, 'T': 0.25,
                     'G': 0.25, 'C': 'a'}
        fragmentation(fasta_dict, seq_counts, nuc_probs)
