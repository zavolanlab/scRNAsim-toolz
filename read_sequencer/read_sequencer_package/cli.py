import argparse
from modules import read_sequencer as rs
import logging

parser = argparse.ArgumentParser(prog='read_sequencer',
                                 description='Simulates sequencing of DNA sequences specified by an FASTA file.')
parser.add_argument('--input_file_path',
                    help='path to FASTA file')
parser.add_argument('--output_file_path',
                    help='path to FASTA file')
parser.add_argument('--read_length',
                    help='read length for sequencing',
                    type=int)

args = parser.parse_args()

def main():
    LOG.info("Program started.")
    read_sequencer = rs()
    read_sequencer.read_fasta(args.input_file_path)
    read_sequencer.run_sequencing(args.read_length)
    read_sequencer.write_fasta(args.output_file_path)
    LOG.info("Program finished.")

if __name__ == '__main__':
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO,
    )
    LOG = logging.getLogger(__name__)
    main()
