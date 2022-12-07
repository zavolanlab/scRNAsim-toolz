import pytest

from readsequencer.read_sequencer import ReadSequencer

sequencer = ReadSequencer()


def test_chunksize():
    assert sequencer.chunk_size == 10000



def test_run_Input():
    assert ReadSequencer(
            fasta="./tests/fasta_testfile/50_seqs_50_1000_bp.fasta",
            output="./tests/fasta_testfile/",
            read_length=1000,
            chunk_size=10000,
        )


def test_run_Random():
    assert ReadSequencer(
            output="./tests/fasta_testfile/",
            read_length=1000,
            chunk_size=10000,
        )