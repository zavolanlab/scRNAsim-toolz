import argparse

parser = argparse.ArgumentParser(prog= 'read_sequencer', description='Simulates Sequenceing of a FASTA file.')
parser.add_argument('--file_path',
                    help='path to FASTA file', action='store_const')
parser.add_argument('--read_length',
                    help='read length for sequencing', action='store_const')


args = parser.parse_args()
print(args.file_path, args.read_length)
