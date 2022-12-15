"""Test utils.py functions."""
import argparse
import pytest

from term_frag_sel.utils import check_positive  # type: ignore


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
