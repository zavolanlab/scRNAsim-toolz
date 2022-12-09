"""Receive command line arguments, fragment sequences, and output fragments."""
import argparse
import logging
from pathlib import Path
from Bio import SeqIO  # type: ignore
import numpy as np
import pandas as pd  # type: ignore


from term_frag_sel.fragmentation import fragmentation
from term_frag_sel.utils import check_positive, check_prob


def main(args: argparse.Namespace):
    """Use CLI arguments to fragment sequences and output text file \
    with selected terminal fragments.

    Args:
        args (parser): list of arguments from CLI.
    """
    # Create or wipe output file
    with open(args.output, "w", encoding="utf-8") as _:
        pass

    logger.info("Checking validity of files...")
    fasta, seq_counts = file_validation(args.fasta, args.counts, args.sep)

    logger.info("Fragmentation of %s...", args.fasta)
    fasta_parse = {}
    for record in fasta:
        fasta_parse[record.id] = record.seq
    splits = np.arange(0, len(list(fasta_parse))+args.size, args.size)

    for i, split in enumerate(splits):
        fasta_dict = fasta_parse[split:splits[i+1]]
        term_frags = fragmentation(fasta_dict, seq_counts,
                                   args.mean, args.std,
                                   args.A_prob, args.T_prob,
                                   args.G_prob, args.C_prob)

        logger.info("Writing batch %s sequences to %s...", i, args.output)
        with open(args.output, 'a', encoding="utf-8") as out_file:
            for line in term_frags:
                out_file.write(f"{line}\n")


def file_validation(fasta_file: str,
                    counts_file: str,
                    sep: str) -> tuple[dict, pd.DataFrame]:
    """Validate input files exist and are the correct format.

    Args:
        fasta_file (str): Input FASTA file path
        counts_file (str): CSV or TSV counts file path
        sep (str): Separator for counts file.

    Returns:
        tuple: fasta and sequence counts variables
    """
    with open(fasta_file, "r", encoding="utf-8") as handle:
        fasta = SeqIO.parse(handle, "fasta")
    if not any(fasta):
        raise ValueError("Input FASTA file is either empty or \
            incorrect file type.")

    count_path = Path(counts_file)
    if not count_path.is_file():
        logger.exception("Input counts file does not exist or isn't a file.")
    else:
        if sep == ",":
            seq_counts = pd.read_csv(counts_file, names=["seqID", "count"])
        else:
            seq_counts = pd.read_table(counts_file, names=["seqID", "count"])

    return fasta, seq_counts


def parse_arguments() -> argparse.Namespace:
    """Request parameters from user on CLI.

    Returns:
        argparse.Namespace: object of arguments from CLI.
    """
    parser = argparse.ArgumentParser(description="""Takes as input FASTA file
                                     of cDNA sequences, a CSV/TSV with sequence
                                     counts, and mean and std. dev. of fragment
                                     lengths and 4 nucleotide probabilities
                                     for the cuts. Outputs most terminal
                                     fragment (within desired length range)
                                     for each sequence.""")

    parser.add_argument('--fasta', required=True,
                        help="Path to FASTA file with cDNA sequences")
    parser.add_argument('--counts', required=True,
                        help="Path to CSV/TSV file with sequence counts")
    parser.add_argument('-o', '--output', required=True,
                        help="output file path")
    parser.add_argument('--mean', required=False, default=300,
                        type=check_positive,
                        help="Mean fragment length (default: 10)")
    parser.add_argument('--std', required=False, default=60,
                        type=check_positive,
                        help="Standard deviation fragment length \
                            (defafult: 1)")
    parser.add_argument('-a', '--A_prob', required=False, default=0.22,
                        type=check_prob,
                        help="Probability cut happens after nucleotide A")
    parser.add_argument('-t', '--T_prob', required=False, default=0.25,
                        type=check_prob,
                        help="Probability cut happens after nucleotide T")
    parser.add_argument('-g', '--G_prob', required=False, default=0.25,
                        type=check_prob,
                        help="Probability cut happens after nucleotide G")
    parser.add_argument('-c', '--C_prob', required=False, default=0.28,
                        type=check_prob,
                        help="Probability cut happens after nucleotide C")
    parser.add_argument('-s', '--size', required=False, default=10000,
                        type=check_positive,
                        help="Chunk size for batch processing")
    parser.add_argument('--sep', required=False, default=",",
                        type=check_positive,
                        help="Sequence counts file separator.")
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s \
            (module "%(module)s")',
        level=logging.INFO,
    )
    logger = logging.getLogger(__name__)

    arguments = parse_arguments()
    main(arguments)
