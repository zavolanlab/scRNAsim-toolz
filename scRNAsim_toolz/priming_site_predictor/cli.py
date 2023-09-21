"""Receive command line arguments."""
import argparse
import logging
from .psp import CreatePrimer
from .psp import PrimingSitePredictor
from scRNAsim_toolz.version import __version__

LOG = logging.getLogger(__name__)


def main():
    """Use CLI arguments to predict priming sites."""
    args = parse_args()
    setup_logging()
    LOG.info("Starting Priming Site Predictor...")

    if args.primer_sequence is None:
        primer = CreatePrimer().create_fasta()
        args.primer_sequence = primer

    predicted_sites = PrimingSitePredictor(
        args.fasta_file, args.primer_sequence, args.energy_cutoff,
        args.riblast_output, args.output_filename
        )

    predicted_sites.generate_gtf()

    LOG.info("Priming Site Predictor finished.")


def setup_logging() -> None:
    """Configure logging.

    Args:
        verbosity: Level of logging verbosity.
    """
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s '
               '(module "%(module)s")',
        level=logging.INFO,
    )


def parse_args():
    """Parse arguments for CLI."""
    parser = argparse.ArgumentParser(
        description="Compute potential priming sites using RIBlast."
        )
    parser.add_argument(
        "-f", "--fasta-file",
        help="Fasta-formatted file of transcript sequences"
        )
    parser.add_argument(
        "-p", "--primer-sequence", default=None,
        help="Primer sequence"
        )
    parser.add_argument(
        "-e", "--energy-cutoff", type=float,
        help="Energy cutoff for interactions"
        )
    parser.add_argument(
        "-r", "--riblast-output",
        help="Path to RIBlast output file"
        )
    parser.add_argument(
        "-o", "--output-filename",
        help="Path where the output gtf should be written"
        )
    parser.add_argument(
        '-v', '--version', action='version',
        version=f'scRNAsim version: {__version__}'
    )
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()
