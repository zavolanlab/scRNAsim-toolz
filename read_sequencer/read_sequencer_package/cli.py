import argparse
from modules import read_sequencer as rs

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
    read_sequencer = rs()
    read_sequencer.read_fasta(args.input_file_path)
    read_sequencer.run_sequencing(args.read_length)
    read_sequencer.write_fasta(args.output_file_path)

if __name__ == '__main__':
    main()
