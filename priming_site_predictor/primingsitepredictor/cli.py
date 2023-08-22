"""Receive command line arguments."""
import argparse
import logging
from primingsitepredictor.prime_site_predictor import CreatePrimer
from primingsitepredictor.prime_site_predictor import PrimingSitePredictor

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
        format='[%(asctime)s: %(levelname)s] %(message)s \
            (module "%(module)s")',
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
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()


# def main():
#     """Execute generate_RIBlast."""
#     generate_riblast_input()


# def generate_riblast_input():
#     """Create a list of the filenames for the RIBlast."""
#     my_primer = CreatePrimer()
#     my_primer.create_fasta()
#     primer_filename = my_primer.name + ".fasta"
#     transcripts_filename = "transcripts.fasta"

#     return [primer_filename, transcripts_filename]


# def create_gtf():
#     """Create gtf."""
#     gtf_file = PostProcessRIBlast().output
#     print(gtf_file)


# def create_parser():
#     """Create a parser."""
#     parser = argparse.ArgumentParser(
#         prog='PrimingSitePredictor',
#         description='Takes a cutoff energy and the predicts \
#             location of priming sites of transcripts',
#         epilog='To predict or not to predict')
#     parser.add_argument('--float', type=float, required=True,
#                         help='An energy-cutoff float number')
#     parsed_args = parser.parse_args()
#     energy_cutoff = parsed_args.float
#     return energy_cutoff


# def letsgo():
#     """Create a parser and print the energycutoff."""
#     energy_cutoff = create_parser()
#     print(f"Your energy cutoff is {energy_cutoff}")
