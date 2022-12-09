"""Test utils.py functions."""
<<<<<<< HEAD
import pytest

from 
=======
import argparse
import pytest

from term_frag_sel.utils import check_positive, check_prob  # type: ignore


def test_positive():
    """Test check_positive function."""
    assert check_positive("1") == 1
    assert check_positive("100") == 100
    assert check_positive("2.2") == 2
    assert check_positive("2.6") == 3
    with pytest.raises(argparse.ArgumentTypeError):
        check_positive("0")
    with pytest.raises(argparse.ArgumentTypeError):
        check_positive("-1")
    with pytest.raises(argparse.ArgumentTypeError):
        check_positive("string")
    with pytest.raises(argparse.ArgumentTypeError):
        check_positive("")


def test_prob():
    """Test check_prob function."""
    assert check_prob("0.1") == 0.1
    assert check_prob("1") == 1.0
    with pytest.raises(argparse.ArgumentTypeError):
        check_prob("0")
    with pytest.raises(argparse.ArgumentTypeError):
        check_prob("10")
    with pytest.raises(argparse.ArgumentTypeError):
        check_prob("-1")
    with pytest.raises(ValueError):
        check_prob("string")
    with pytest.raises(ValueError):
        check_prob("")
>>>>>>> hugo_new
