"""Command line interface."""
import argparse
import logging
from pathlib import Path

from tsg.main import sample_transcripts


def setup_logging(loglevel: str) -> None:
    """Set up logging. Loglevel can be one of \
        ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"].

    Args:
        loglevel: Level of log output.

    Raises:
        ValueError: If string that is not a log level is passed, raise error.

    Returns:
        None
    """
    # set default log level
    numeric_level = getattr(logging, "INFO")

    if loglevel:
        try:
            numeric_level = getattr(logging, loglevel.upper())
        except AttributeError as err:
            print(f"Unexpected {err=}, {type(err)=}")
            raise

    logging.basicConfig(
        format=('[%(asctime)s: %(levelname)s] '
                '%(message)s (module "%(module)s")'),
        level=numeric_level,
    )


def build_arg_parser() -> argparse.Namespace:
    """Build the argument parser.

    Args:
        1) path to the csv-file with the number of transcripts
        2) path to the gtf-file with the annotations for each transcript
        3) a value for the probability of intron inclusion
        4) a log message

    Raises:
        None

    Returns:
        arguments for parser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "transcripts",
        type=str,
        help="Path to csv file with number of transcripts (ID,Count).",
    )
    parser.add_argument(
        "annotation",
        type=str,
        help="Path to gtf-file with exon annotation."
    )
    parser.add_argument(
        "-p",
        "--prob-inclusion",
        type=float,
        default=0.05,
        help="Probability of intron inclusion.",
    )
    parser.add_argument(
        "--log",
        type=str,
        default="INFO",
        help='Level of logging. Can be one of \
            ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]',
    )

    args = parser.parse_args()

    assert args.prob_inclusion >= 0
    assert args.prob_inclusion <= 1

    return args


def output_filename(filename: str) -> str:
    """Generate output filename for given input filename.

    Args:
        filename: Input filename

    Raises:
        NotImplementedError: Only accept filetypes .csv, .tsv and .gtf.
        FileExistsError: If the output file exists, raise error.

    Returns:
        str: Output filename
    """
    filepath = Path(filename)
    if filepath.suffix in (".csv", ".tsv"):
        outfile = "generated_" + filepath.stem + ".csv"
    elif filepath.suffix == ".gtf":
        outfile = "generated_" + filepath.name
    else:
        raise NotImplementedError()

    if Path(outfile).exists():
        raise FileExistsError(f"The output file {outfile} already exists.")

    return outfile


def app():
    """Get the args, sets up the logging \
        and starts the programm with the provided parameters.

    Args:
        1) path to the csv-file with the number of transcripts
        2) path to the gtf-file with the annotations for each transcript
        3) a value for the probability of intron inclusion
        4) a log message

    Raises:
        None

    Returns:
        None
    """
    args = build_arg_parser()

    setup_logging(args.log)
    sample_transcripts(
        args.transcripts,
        args.annotation,
        args.prob_inclusion,
        output_filename(args.transcripts),
        output_filename(args.annotation),
    )
