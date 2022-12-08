"""Utility functions for command line arguments."""
import argparse
import os.path


# found on shorturl.at/vzAX4
def extant_file(x):
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    elif not x.endswith((".fasta", ".fa", ".csv")):
        raise argparse.ArgumentTypeError("""{0} is not the correct
                                         file format""".format(x))
    return x


def check_positive(value):
    """Check input value is a positive integer.

    Args:
        value (string): command line parameter

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
    else:
        return ivalue


def check_prob(value):
    """
    Check probability value is within ]0,1] range.

    Args:
        value (string): command line parameter

    Raises:
        argparse.ArgumentTypeError: received a value outside valid range

    Returns:
        float: float version of input value
    """
    pvalue = float(value)
    if pvalue <= 0 or pvalue > 1:
        raise argparse.ArgumentTypeError("""%s is an invalid positive int
                                         value""" % value)
    return pvalue
