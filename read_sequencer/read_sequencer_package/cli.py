import argparse
from modules import ReadSequencer
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
parser.add_argument('--random', action='store_true', default=False,
                    help='generate random sequences')
parser.add_argument('--n_random', default=100, type=int, help='n random sequences')
parser.add_argument('--mean_random', default=50, type=int, help='mean random sequences')
parser.add_argument('--sd_random', default=25, type=int, help='standard deviation random sequences')

args = parser.parse_args()

def main():
    LOG.info("Read sequencer started.")
    read_sequencer = ReadSequencer()
    if args.random:
        read_sequencer.add_random_sequences(n=args.n_random, mean=args.mean_random, sd=args.sd_random)
    else:
        read_sequencer.read_fasta(args.input_file_path)
    read_sequencer.run_sequencing(args.read_length)
    read_sequencer.write_fasta(args.output_file_path)
    LOG.info("Read sequencer finished.")

if __name__ == '__main__':
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO)
    LOG = logging.getLogger(__name__)
    main()
