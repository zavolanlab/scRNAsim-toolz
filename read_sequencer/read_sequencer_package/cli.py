import argparse
from modules import run_read_sequencer

parser = argparse.ArgumentParser(prog='read_sequencer',
                                 description='Simulates sequencing of a DNA sequences specified by an FASTA file.')
parser.add_argument('--input_file_path',
                    help='path to FASTA file')
parser.add_argument('--output_file_path',
                    help='path to FASTA file')
parser.add_argument('--read_length',
                    help='read length for sequencing',
                    type=int)

args = parser.parse_args()


def main():
    run_read_sequencer(args.input_file_path, args.read_length, args.output_file_path)


if __name__ == '__main__':
    main()
