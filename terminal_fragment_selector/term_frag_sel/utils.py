"""Utility functions for CLI arguments."""
import argparse


def check_positive(value: str) -> int:
    """Check input value is a positive integer.

    Args:
        value (str): command line parameter

    Raises:
        argparse.ArgumentTypeError: received a negative integer
        argparse.ArgumentTypeError: received a non-integer value

    Returns:
        integer: integer version of input value
    """
    try:
        ivalue = round(float(value))
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(f"""Expected positive integer,
                                             got negative integer: {value}""")
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"""Expected positive integer,
                                         got: {value}""") from exc
    return ivalue
