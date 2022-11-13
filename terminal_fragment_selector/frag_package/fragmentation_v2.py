import argparse
import os.path
import re

import numpy as np


def fragmentation(fasta, seq_counts, mean_length, std):
    nucs = ['A','T','G','C']
    mononuc_freqs = [0.22, 0.25, 0.23, 0.30]
    term_frags = [] 
    for seq, counts in dna_seq.items():
        for _ in range(counts): 
            n_cuts = int(len(seq)/mean_length)
            
            # non-uniformly random DNA fragmentation implementation based on https://www.nature.com/articles/srep04532#Sec1
            # assume fragmentation by sonication for NGS workflow
            cuts = []
            cut_nucs = np.random.choice(nucs, n_cuts, p=mononuc_freqs) 
            for nuc in cut_nucs:
                nuc_pos = [x.start() for x in re.finditer(nuc, seq)]
                pos = np.random.choice(nuc_pos)
                while pos in cuts:
                    pos = np.random.choice(nuc_pos)
                cuts.append(pos)

            cuts.sort() 
            cuts.insert(0,0)
            term_frag = ""
            for i, val in enumerate(cuts):
                if i == len(cuts)-1:
                    fragment = seq[val+1:cuts[-1]]
                else:
                    fragment = seq[val:cuts[i+1]]
                if mean_length-std <= len(fragment) <= mean_length+std:
                    term_frag = fragment
            if term_frag == "":
                continue
            else:
                term_frags.append(term_frag)
    return term_frags

def main(args):   
    fasta, seq_counts, mean_length, std = args
    dna_seq = {
    "ATAACATGTGGATGGCCAGTGGTCGGTTGTTACACGCCTACCGCGATGCTGAATGACCCGGACTAGAGTGGCGAAATTTATGGCGTGTGACCCGTTATGC": 100,
    "TCCATTTCGGTCAGTGGGTCATTGCTAGTAGTCGATTGCATTGCCATTCTCCGAGTGATTTAGCGTGACAGCCGCAGGGAACCCATAAAATGCAATCGTA": 100}

    term_frags = fragmentation(fasta, seq_counts, mean_length, std)
    with open('terminal_frags.txt', 'w') as f:
        for line in term_frags:
            f.write(line)
            f.write('\n')

# found on https://stackoverflow.com/questions/11540854/file-as-command-line-argument-for-argparse-error-message-if-argument-is-not-va
def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Takes as input FASTA file of cDNA sequences, a CSV with sequence counts, and mean and std. dev. of fragment lengths. Outputs most terminal fragment (within desired length range) for each sequence.")
    
    parser.add_argument('--fasta', required=True, type=extant_file, help="FASTA file with cDNA sequences")
    parser.add_argument('--counts', required=True, type=extant_file, help="CSV file with sequence counts")
    parser.add_argument('--mean', required = False, default = 10, type = int, help="Mean fragment length (default: 10)")
    parser.add_argument('--std', required = False, default = 1, type = int, help="Standard deviation fragment length (defafult: 1)")
    args = parser.parse_args()
    
    return args.fasta, args.counts, args.mean, args.std


if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)
