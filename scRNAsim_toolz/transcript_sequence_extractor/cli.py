"""CLI to be run on output fasta file from bedtools getfasta."""
import argparse
import logging
from sequence_extractor.pre_bedtools import pre_bedtools_mode
from sequence_extractor.exon_concatenation import exon_concatenation
from sequence_extractor.poly_a import poly_a_addition_to_fasta_list

LOG = logging.getLogger(__name__)


def main():
    """Use CLI arguments to extract sequences.

    Runs on the output from bedtools and concatenates the exons together
    and adds a polyA tail and outputs a fasta file.

    Args:
        None: this will run on its own by taking the information from argparse

    Returns:
        A fasta file with a single entry for each transcript ID with
        polyA tail being added onto the sequence at 3'end
    """
    args = parse_args()
    setup_logging()

    if args.mode == "pre_bedtools":
        pre_bedtools_mode(args)
    elif args.mode == "post_bedtools":
        post_bedtools_mode(args)
    else:
        LOG.error(
            "Invalid mode specified."
            "Please choose 'pre_bedtools' or 'post_bedtools'.")


def setup_logging() -> None:
    """Configure logging."""
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s '
               '(module "%(module)s")',
        level=logging.INFO,
    )


def parse_args():
    """Parse arguments for CLI."""
    parser = argparse.ArgumentParser(
        description="extracts transcript sequences from genome sequence and"
                    "ouputs transcripts with PolyA tail added to them",
        )
    parser.add_argument(
        "--mode",
        choices=["pre_bedtools", "post_bedtools"],
        required=True,
        help="Select the mode of operation"
             "('pre_bedtools' or 'post_bedtools')."
    )
    parser.add_argument(
        "-i", "--input-fasta-file",
        dest="input_fasta_file",
        help="Fasta-formatted file obtained from bedtools"
    )
    parser.add_argument(
        "-o", "--output-file-name",
        dest="output_file_name",
        help="Name of the output fasta file"
    )
    parser.add_argument(
        "-p", "--polyA-length",
        dest="poly_a_length",
        type=int,
        help="Length of the polyA tail to be added (def: 250)",
        default=250
    )
    parser.add_argument(
        "--input-gtf-file",
        dest="input_gtf_file",
        help="Ordered and processed gtf file for 'pre_bedtools' mode.")
    parser.add_argument(
        "--output-bed-file",
        dest="output_bed_file",
        help="Bed file with only exons with strandedness"
             "taken into account for 'pre_bedtools' mode.")

    args = parser.parse_args()
    return args


def post_bedtools_mode(args):
    """Execute the 'post_bedtools' mode."""
    LOG.info("Starting 'post_bedtools' mode...")

    fasta_list = exon_concatenation(args.input_fasta_file)
    final_list = poly_a_addition_to_fasta_list(fasta_list, args.poly_a_length)

    if args.output_file_name is None:
        args.output_file_name = "default_output.fasta"

    with open(args.output_file_name, "w", encoding="utf-8") as fasta_out:
        fasta_out.write("\n".join(f"{x[0]}\n{x[1]}" for x in final_list))
    LOG.info("Transcript sequence extractor finished in 'post_bedtools' mode.")


if __name__ == '__main__':
    main()
