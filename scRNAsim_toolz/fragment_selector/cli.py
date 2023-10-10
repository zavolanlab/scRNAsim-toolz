"""Receive command line arguments, fragment sequences, and output fragments."""
import argparse
import logging
from pathlib import Path
from Bio import SeqIO  # type: ignore
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from scRNAsim_toolz.version import __version__

from .fragmentation import fragmentation
from .utils import check_positive

logging.basicConfig(
    format='[%(asctime)s: %(levelname)s] %(message)s \
        (module "%(module)s")',
    level=logging.INFO,
)
logger = logging.getLogger("main")


def main():
    """Use CLI arguments to fragment sequences and output text file \
    with selected terminal fragments.

    Args:
        args (parser): list of arguments from CLI.
    """
    args = parse_arguments()

    if not isinstance(args, argparse.Namespace):
        raise TypeError("Input should be argparse.Namespace")

    # Create or wipe output file
    with open(args.output, "w", encoding="utf-8") as _:
        pass

    logger.info("Checking validity of files...")
    fasta, seq_counts = file_validation(args.fasta, args.counts, args.sep)

    logger.info("Fragmentation of %s...", args.fasta)
    splits = np.arange(0, len(list(fasta))+args.size, args.size)

    for i in range(len(splits) - 1):
        split = splits[i]
        keys = list(fasta.keys())[split:splits[i+1]]
        fasta_dict = {key: fasta[key] for key in keys}
        term_frags = fragmentation(fasta_dict, seq_counts,
                                   args.mean, args.std)

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
        tuple: fasta dict and sequence counts pd.DataFrame
    """
    fasta_dict = {}
    with open(fasta_file, "r", encoding="utf-8") as handle:
        fasta_sequences = SeqIO.parse(handle, "fasta")

        if not any(fasta_sequences):
            raise ValueError("Input FASTA file is either empty or \
                incorrect file type.")

        for record in fasta_sequences:
            fasta_dict[record.id] = str(record.seq).upper()

    count_path = Path(counts_file)
    if not count_path.is_file():
        raise FileNotFoundError("Input counts file does not exist or \
            isn't a file.")

    if sep == ",":
        seq_counts = pd.read_csv(counts_file, names=["seqID", "count"])
        seq_counts = seq_counts.astype({"seqID": str})
    else:
        seq_counts = pd.read_table(counts_file, names=["seqID", "count"])
        seq_counts = seq_counts.astype({"seqID": str})

    return fasta_dict, seq_counts


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
                        help="Mean fragment length (default: 300)")
    parser.add_argument('--std', required=False, default=60,
                        type=check_positive,
                        help="Standard deviation fragment length \
                            (defafult: 60)")
    parser.add_argument('-s', '--size', required=False, default=10000,
                        type=check_positive,
                        help="Chunk size for batch processing")
    parser.add_argument('--sep', required=False, default=",",
                        help="Sequence counts file separator.")
    parser.add_argument(
        '-v', '--version', action='version',
        version=f'scRNAsim version: {__version__}'
    )
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s \
            (module "%(module)s")',
        level=logging.INFO,
    )
    logger = logging.getLogger(__name__)
    main()
