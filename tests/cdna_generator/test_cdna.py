"""Tests for cDNA functions."""
import pytest  # type: ignore
from ...scRNAsim_toolz.cdna_generator.cdna import (  # type: ignore
    complement, seq_complement
)


@pytest.mark.parametrize(
    "test_input,expected",
    [("A", "T")]
)
# we need to pass the lists to the test function...
def test_complement_param(test_input, expected):
    """Test complement() function."""
    assert complement(test_input) == expected


@pytest.mark.parametrize(
    "test_input,expected",
    [("AA", "TT")]
)
# we need to pass the lists to the test function...
def test_seq_complement_param(test_input, expected):
    """Test seq_complement() function."""
    assert seq_complement(test_input) == expected


# we can do the same for the tests that raise an error:
@pytest.mark.parametrize(
    "test_input,expected",
    [(1, ValueError)]
)
def test_complement_param_failing(test_input, expected):
    """Test complement() fail function."""
    with pytest.raises(expected):
        complement(test_input)


@pytest.mark.parametrize(
    "test_input,expected",
    [("11", ValueError)]
)
def test_seq_complement_param_failing(test_input, expected):
    """Test seq_complement() fail function."""
    with pytest.raises(expected):
        seq_complement(test_input)
