import pytest
import os
import glob
from readsequencer.read_sequencer import ReadSequencer

def test_init_default():
    sequencer = ReadSequencer()
    assert sequencer.fasta is None
    assert sequencer.read_length == 150
    assert sequencer.output is None
    assert sequencer.chunk_size == 10000
    assert sequencer.bases == ("A", "T", "C", "G")


def test_run_random():
    sequencer = ReadSequencer(
        output="./tests/fasta_testfile/results.fasta")
    sequencer.define_random_sequences(n_seq=100)
    assert sequencer.output == "./tests/fasta_testfile/results.fasta"
    assert sequencer.read_length == 150
    assert sequencer.chunk_size == 10000
    assert sequencer.fasta is None
    sequencer.run_sequencing()
    os.remove("./tests/fasta_testfile/results.fasta")

def test_run_random_chunks():
    # setup class
    sequencer = ReadSequencer(
        output="./tests/fasta_testfile/results.fasta",
        read_length=150,
        chunk_size=10)
    sequencer.define_random_sequences(n_seq=50)
    # run sequencing
    sequencer.run_sequencing()
    # check results
    assert sequencer.output == "./tests/fasta_testfile/results.fasta"
    assert sequencer.read_length == 150
    assert sequencer.n_sequences == 50
    # clean up
    result_files = glob.glob("./tests/fasta_testfile/results*")
    assert len(result_files) == 5
    for file in result_files:
        os.remove(file)


def test_run_sequencing():
    sequencer = ReadSequencer(
        fasta="./tests/fasta_testfile/50_seqs_50_1000_bp.fasta",
        output="./tests/fasta_testfile/results.fasta",
        read_length=50,
        chunk_size=10000)
    sequencer.get_n_sequences()
    sequencer.run_sequencing()
    assert sequencer.output == "./tests/fasta_testfile/results.fasta"
    assert sequencer.read_length == 50
    assert sequencer.n_sequences == 50
    result_file = glob.glob("./tests/fasta_testfile/results*")
    assert len(result_file) == 1
    for file in result_file:
        os.remove(file)

def test_run_sequencing_chunks():
    # setup class
    sequencer = ReadSequencer(
        fasta="./tests/fasta_testfile/50_seqs_50_1000_bp.fasta",
        output="./tests/fasta_testfile/results.fasta",
        read_length=150,
        chunk_size=10)
    # run sequencing
    sequencer.get_n_sequences()
    sequencer.run_sequencing()
    # check results
    assert sequencer.output == "./tests/fasta_testfile/results.fasta"
    assert sequencer.read_length == 150
    assert sequencer.n_sequences == 50
    # clean up
    result_files = glob.glob("./tests/fasta_testfile/results*")
    assert len(result_files) == 5
    for file in result_files:
        os.remove(file)



