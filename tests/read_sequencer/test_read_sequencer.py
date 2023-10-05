"""Test read_sequencer.py."""
import os
import glob
from pathlib import Path
from scRNAsim_toolz.read_sequencer.read_sequencer import (
    ReadSequencer
)

TEST_FILES_DIR = Path(__file__).resolve().parent / "fasta_testfile"


def test_init_default():
    """Test default initation."""
    sequencer = ReadSequencer()
    assert sequencer.fasta is None
    assert sequencer.read_length == 150
    assert sequencer.output is None
    assert sequencer.chunk_size == 10000
    assert sequencer.bases == ("A", "T", "C", "G")


def test_run_random():
    """Test random run."""
    sequencer = ReadSequencer(
        output=str(TEST_FILES_DIR / "results.fasta")
    )
    sequencer.define_random_sequences(n_seq=100)
    assert sequencer.output == str(TEST_FILES_DIR / "results.fasta")
    assert sequencer.read_length == 150
    assert sequencer.chunk_size == 10000
    assert sequencer.fasta is None
    sequencer.run_sequencing()
    os.remove(TEST_FILES_DIR / "results.fasta")


def test_run_random_chunks():
    """Test random run chunks."""
    # setup class
    sequencer = ReadSequencer(
        output=str(TEST_FILES_DIR / "results.fasta"),
        read_length=150,
        chunk_size=10)
    sequencer.define_random_sequences(n_seq=50)
    # run sequencing
    sequencer.run_sequencing()
    # check results
    assert sequencer.output == str(TEST_FILES_DIR / "results.fasta")
    assert sequencer.read_length == 150
    assert sequencer.n_sequences == 50
    # clean up
    result_files = list(TEST_FILES_DIR.glob("results*"))
    assert len(result_files) == 5
    for file in result_files:
        os.remove(file)


def test_run_sequencing():
    """Test sequencing run."""
    sequencer = ReadSequencer(
        fasta=str(TEST_FILES_DIR / "50_seqs_50_1000_bp.fasta"),
        output=str(TEST_FILES_DIR / "results.fasta"),
        read_length=50,
        chunk_size=10000)
    sequencer.get_n_sequences()
    sequencer.run_sequencing()
    assert sequencer.output == str(TEST_FILES_DIR / "results.fasta")
    assert sequencer.read_length == 50
    assert sequencer.n_sequences == 50
    result_file = list(TEST_FILES_DIR.glob("results*"))
    assert len(result_file) == 1
    for file in result_file:
        os.remove(file)


def test_run_sequencing_chunks():
    """Test run sequencing chunks."""
    # setup class
    sequencer = ReadSequencer(
        fasta=str(TEST_FILES_DIR / "50_seqs_50_1000_bp.fasta"),
        output=str(TEST_FILES_DIR / "results.fasta"),
        read_length=150,
        chunk_size=10)
    # run sequencing
    sequencer.get_n_sequences()
    sequencer.run_sequencing()
    # check results
    assert sequencer.output == str(TEST_FILES_DIR / "results.fasta")
    assert sequencer.read_length == 150
    assert sequencer.n_sequences == 50
    # clean up
    result_files = list(TEST_FILES_DIR.glob("results*"))
    assert len(result_files) == 5
    for file in result_files:
        os.remove(file)
