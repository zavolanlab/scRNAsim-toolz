"""Fragment sequences."""
import re

import numpy as np
import pandas as pd


def fasta_process(fasta_file):
    """
    Pre-process FASTA file.

    Args:
        fasta_file (fasta): FASTA file with cDNA sequences

    Returns:
        dict: Dictionary of gene sequence IDs and their sequence
    """
    with open(fasta_file, "r") as f:
        lines = f.readlines()

        # Tanya, try \\S instead of \S and see if that works
        ident_pattern = re.compile('>(\S+)')
        seq_pattern = re.compile('^(\S+)$')

        genes = {}
        for line in lines:
            if ident_pattern.search(line):
                seq_id = (ident_pattern.search(line)).group(1)
            elif seq_id in genes.keys():
                genes[seq_id] += (seq_pattern.search(line)).group(1)
            else:
                genes[seq_id] = (seq_pattern.search(line)).group(1)
    return genes


def fragmentation(fasta_file, counts_file, mean_length, std,
                  a_prob, t_prob, g_prob, c_prob):
    """
    Fragment cDNA sequences and select terminal fragment.

    Args:
        fasta_file (fasta): FASTA file with cDNA sequences
        counts_file (text): CSV or TSV file woth sequence counts
        mean_length (int): mean length of desired fragments
        std (int): standard deviation of desired fragment lengths
        a_prob (float): probability of nucleotide A
        t_prob (float): probability of nucleotide T
        g_prob (float): probability of nucleotide G
        c_prob (float): probability of nucleotide C

    Returns:
        list: list of selected terminal fragments
    """
    fasta = fasta_process(fasta_file)
    seq_counts = pd.read_csv(counts_file,
                             names=["seqID", "count"])

    # calculated using https://www.nature.com/articles/srep04532#MOESM1
    nuc_probs = {'A': a_prob, 'T': t_prob, 'G': g_prob, 'C': c_prob}

    term_frags = []
    for seq_id, seq in fasta.items():
        counts = seq_counts[seq_counts["seqID"] == seq_id]["count"]
        for _ in range(counts):
            n_cuts = int(len(seq)/mean_length)

            # non-uniformly random DNA fragmentation implementation based on
            # https://www.nature.com/articles/srep04532#Sec1
            # assume fragmentation by sonication for NGS workflow
            cuts = []
            cut_nucs = np.random.choice(list(nuc_probs.keys()),
                                        n_cuts, p=list(nuc_probs.values()))
            for nuc in cut_nucs:
                nuc_pos = [x.start() for x in re.finditer(nuc, seq)]
                pos = np.random.choice(nuc_pos)
                while pos in cuts:
                    pos = np.random.choice(nuc_pos)
                cuts.append(pos)

            cuts.sort()
            cuts.insert(0, 0)
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
