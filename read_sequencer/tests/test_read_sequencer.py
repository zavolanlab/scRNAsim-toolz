import pytest

from readsequencer.read_sequencer import ReadSequencer

sequencer = ReadSequencer()


def test_chunksize():
    assert sequencer.chunk_size == 10000
