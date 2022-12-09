import pytest

from readsequencer.read_sequencer import ReadSequencer


sequencer = ReadSequencer()


def test_chunksize():
    assert sequencer.chunk_size == 10000


def test_run_Input():
    assert sequencer.fasta == None
    assert sequencer.read_length == 150
    assert sequencer.output == None
    assert sequencer.chunk_size == 10000


def test_run_Random():
    assert ReadSequencer(
        output="./tests/fasta_testfile/results.fasta"
    ).output == "./tests/fasta_testfile/results.fasta"
    assert ReadSequencer(read_length=1000).read_length == 1000
    assert ReadSequencer(chunk_size=10000).chunk_size == 10000
    assert ReadSequencer(
        output="./tests/fasta_testfile/results.fasta",
        read_length=1000,
        chunk_size=10000).fasta == None
