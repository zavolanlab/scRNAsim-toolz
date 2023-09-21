"""Receive command line arguments."""
import argparse
import logging
from scRNAsim_toolz.version import __version__

logging.basicConfig(
    format='[%(asctime)s: %(levelname)s] %(message)s \
        (module "%(module)s")',
    level=logging.INFO,
)
LOG = logging.getLogger(__name__)

from .cdna import CDNAGen  # noqa: E402,E501 # pylint:disable=wrong-import-position


def main():
    """Parse sequences for cDNA generator.

    Parses command line arguments for cDNA generation.

    Returns: CDNAGen instance

    """
    parser = argparse.ArgumentParser(
        prog="cDNA generator",
        description="Generate cDNA sequences based on primer probabilities.",
    )
    parser.add_argument(
        "-ifa", "--input_fasta", help="genome fasta file", required=True
    )
    parser.add_argument("-igtf", "--input_gtf", help="gtf file", required=True)
    parser.add_argument(
        "-ofa", "--output_fasta", help="output fasta file", required=True
    )
    parser.add_argument(
        "-icpn",
        "--input_copy_number",
        help="input copy number (csv) file",
        required=True,
    )
    parser.add_argument(
        "-ocsv", "--output_csv", help="output fasta file", required=True
    )
    parser.add_argument(
        '-v', '--version', action='version',
        version=f'scRNAsim version: {__version__}'
    )
    args = parser.parse_args()

    LOG.info("Running cDNA generator...")
    CDNAGen(
        ifasta=args.input_fasta,
        igtf=args.input_gtf,
        icpn=args.input_copy_number,
        ocsv=args.output_csv,
        ofasta=args.output_fasta,
    )


if __name__ == "__main__":
    main()
