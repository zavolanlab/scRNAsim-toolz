"""Fragment sequences."""
import re

import numpy as np
import pandas as pd  # type: ignore


def fragmentation(fasta: dict, seq_counts: pd.DataFrame,
                  mean_length: int = 100, std: int = 10,
                  a_prob: float = 0.22, t_prob: float = 0.25,
                  g_prob: float = 0.25, c_prob: float = 0.28
                  ) -> list:
    """Fragment cDNA sequences and select terminal fragment.

    Args:
        fasta_file (dict): dictionary of {transcript IDs: sequences}
        counts_file (pd.DataFrame): dataframe with sequence counts and IDs
        mean_length (int): mean length of desired fragments
        std (int): standard deviation of desired fragment lengths
        a_prob (float): probability of nucleotide A
        t_prob (float): probability of nucleotide T
        g_prob (float): probability of nucleotide G
        c_prob (float): probability of nucleotide C

    Returns:
        list: list of selected terminal fragments
    """
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
            if term_frag != "":
                term_frags.append(term_frag)

    return term_frags
