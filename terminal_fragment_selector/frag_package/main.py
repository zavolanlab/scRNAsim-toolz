import argparse

from fragmentation_v2 import fragmentation
from utils import check_positive, extant_file


def main(args):   
    fasta, seq_counts, mean_length, std = args

    term_frags = fragmentation(fasta, seq_counts, mean_length, std)
    with open('terminal_frags.txt', 'w') as f:
        for line in term_frags:
            f.write(line)
            f.write('\n')

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Takes as input FASTA file of cDNA sequences, a CSV with sequence counts, and mean and std. dev. of fragment lengths. Outputs most terminal fragment (within desired length range) for each sequence.")

    parser.add_argument('--fasta', required=True, type=extant_file, help="FASTA file with cDNA sequences")
    parser.add_argument('--counts', required=True, type=extant_file, help="CSV file with sequence counts")
    parser.add_argument('--mean', required = False, default = 10, type = check_positive, help="Mean fragment length (default: 10)")
    parser.add_argument('--std', required = False, default = 1, type = check_positive, help="Standard deviation fragment length (defafult: 1)")
    args = parser.parse_args()
    
    return args.fasta, args.counts, args.mean, args.std


if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)