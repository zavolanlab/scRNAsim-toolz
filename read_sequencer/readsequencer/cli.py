import argparse
import logging
from readsequencer.read_sequencer import ReadSequencer

LOG = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        prog="readsequencer",
        description="Simulates sequencing of DNA sequences specified by an FASTA file.",
    )

    parser.add_argument("output", help="path to FASTA file")
    parser.add_argument("-i", "--input", default=None, help="path to FASTA file")
    parser.add_argument(
        "-r", "--read-length", default=100, help="read length for sequencing", type=int
    )
    parser.add_argument(
        "-n",
        "--n_random",
        default=100,
        type=int,
        help="n random sequences. Just used if input fasta file is not specified.",
    )
    parser.add_argument(
        "-s",
        "--chunk-size",
        default=10000,
        type=int,
        help="chunk_size for batch processing",
    )

    args = parser.parse_args()
    LOG.info("Read sequencer started.")
    if args.input is not None:
        read_sequencer = ReadSequencer(
            fasta=args.input,
            output=args.output,
            read_length=args.read_length,
            chunk_size=args.chunk_size,
        )
        read_sequencer.get_n_sequences()
    else:
        read_sequencer = ReadSequencer(
            fasta=args.input,
            output=args.output,
            read_length=args.read_length,
            chunk_size=args.chunk_size,
        )
        read_sequencer.define_random_sequences(n_seq=args.n_random)

    read_sequencer.run_sequencing()

    LOG.info("Read sequencer finished.")


if __name__ == "__main__":
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO,
    )
    LOG = logging.getLogger(__name__)
    main()
