import argparse
import os.path


# found on https://stackoverflow.com/questions/11540854/file-as-command-line-argument-for-argparse-error-message-if-argument-is-not-va
def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    elif not x.endswith((".fasta", ".fa", ".csv")):
        raise argparse.ArgumentTypeError("{0} is not the correct file format".format(x))
    return x

# found on https://stackoverflow.com/questions/14117415/in-python-using-argparse-allow-only-positive-integers
def check_positive(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

def check_prob(value):
    pvalue = float(value)
    if pvalue <= 0 or pvalue>1:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return pvalue 