import argparse
import logging
from pathlib import Path

from tsg.main import sample_transcripts


def setup_logging(loglevel: str = None) -> None:
    # Set up logging
    if loglevel:
        numeric_level = getattr(logging, loglevel.upper())
        if not isinstance(numeric_level, int):
            raise ValueError("Invalid log level: %s" % loglevel)
    else:
        numeric_level = getattr(logging, "INFO")

    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=numeric_level,
    )


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--transcripts", type=str)
    parser.add_argument("--annotation", type=str)
    parser.add_argument("--prob_inclusion", type=float)
    parser.add_argument("--log", type=str)

    return parser


def get_args():
    parser = build_arg_parser()

    args = parser.parse_args()

    return args


def output_filename(filename: str) -> str:
    filepath = Path(filename)
    if filepath.suffix == ".csv" or filepath.suffix == ".tsv":
        outfile = "generated_" + filepath.stem + ".csv"
    elif filepath.suffix == ".gtf":
        outfile = "generated_" + filepath.name
    else:
        raise NotImplementedError()
    return outfile


def app():
    args = get_args()

    setup_logging(args.log)
    sample_transcripts(
        args.transcripts,
        args.annotation,
        args.prob_inclusion,
        output_filename(args.transcripts),
        output_filename(args.annotation),
    )
