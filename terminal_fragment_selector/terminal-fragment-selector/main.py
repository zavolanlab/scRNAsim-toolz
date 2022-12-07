import argparse

from fragmentation_v2 import fragmentation
from utils import check_positive, extant_file, check_prob


def main(args):   
    fasta, seq_counts, mean_length, std, a_prob, t_prob, g_prob, c_prob = args

    term_frags = fragmentation(fasta, seq_counts, mean_length, std, a_prob, t_prob, g_prob, c_prob)
    with open('terminal_frags.txt', 'w') as f:
        for line in term_frags:
            f.write(line)
            f.write('\n')

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Takes as input FASTA file of cDNA sequences, a CSV with sequence counts, and mean and std. dev. of fragment lengths and 4 nucleotide probabilities for the cuts. Outputs most terminal fragment (within desired length range) for each sequence.")

    parser.add_argument('--fasta', required=True, type=extant_file, help="FASTA file with cDNA sequences")
    parser.add_argument('--counts', required=True, type=extant_file, help="CSV file with sequence counts")
    parser.add_argument('--mean', required = False, default = 300, type = check_positive, help="Mean fragment length (default: 10)")
    parser.add_argument('--std', required = False, default = 60, type = check_positive, help="Standard deviation fragment length (defafult: 1)")
    parser.add_argument('--A_prob', required=False, default = 0.22, type=check_prob, help="Probability cut happens after nucleotide A")
    parser.add_argument('--T_prob', required=False, default = 0.25, type=check_prob, help="Probability cut happens after nucleotide T")
    parser.add_argument('--G_prob', required=False, default = 0.25, type=check_prob, help="Probability cut happens after nucleotide G")
    parser.add_argument('--C_prob', required=False, default = 0.28, type=check_prob, help="Probability cut happens after nucleotide C")
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)