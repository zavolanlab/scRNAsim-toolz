# imports
import pytest
from cdna.cdna import compliment, seq_compliment

@pytest.mark.parametrize(
    "test_input,expected",
    [("A", "T")]
)
def test_compliment_param(test_input, expected):  # we need to pass the lists to the test function...
    assert compliment(test_input) == expected

@pytest.mark.parametrize(
    "test_input,expected",
    [("AA", "TT")]
)
def test_seq_compliment_param(test_input, expected):  # we need to pass the lists to the test function...
    assert seq_compliment(test_input) == expected


# we can do the same for the tests that raise an error:
@pytest.mark.parametrize(
    "test_input,expected",
    [(1, ValueError)]
)
def test_compliment_param_failing(test_input, expected):
    with pytest.raises(expected):
        compliment(test_input)

@pytest.mark.parametrize(
    "test_input,expected",
    [("11", ValueError)]
)
def test_compliment_param_failing(test_input, expected):
    with pytest.raises(expected):
        seq_compliment(test_input)
