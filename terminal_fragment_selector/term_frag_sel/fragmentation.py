"""Fragment sequences."""
import re

import numpy as np
import pandas as pd  # type: ignore


def get_cut_number(seq_len: int, mean: float) -> int:
    """Get the number of cuts for a particular sequence.

    Args:
        seq_len (int): length of sequence/number of nucleotides in the sequence
        mean (float): mean fragment length

    Returns:
        int: number of cuts
    """
    if not isinstance(seq_len, int):
        raise ValueError(f"Sequence length must be numeric, \
            not {type(seq_len)}")

    if not isinstance(mean, int):
        raise ValueError(f"Mean must be numeric, not {type(mean)}")

    cuts_distribution = []  # distribution of cut numbers (n_cuts)

    # 1000 iterations should not be too computationally inefficient
    # given the nature of the code
    for _ in range(1000):
        n_cuts = 0
        len_sum = 0
        while True:
            len_sum += np.random.exponential(scale=mean)
            if len_sum < seq_len:
                n_cuts += 1
            else:
                cuts_distribution.append(n_cuts)
                break

    cuts_distribution.sort()
    cut_counts = {x: cuts_distribution.count(x) for x in cuts_distribution}
    cut_probs = [x/1000 for x in cut_counts.values()]

    # randomly ick no. of cut from cut distribution based on probability
    # of cut numbers
    n_cuts = np.random.choice(list(cut_counts.keys()), p=cut_probs)

    return n_cuts


def fragmentation(fasta: dict, seq_counts: pd.DataFrame,
                  nuc_probs: dict,
                  mu_length: int = 100, std: int = 10,
                  ) -> list:
    """Fragment cDNA sequences and select terminal fragment.

    Args:
        fasta_file (dict): dictionary of {transcript IDs: sequences}
        counts_file (pd.DataFrame): dataframe with sequence counts and IDs
        nuc_probs (dict): probability of cut occuring a certain nucleotide. \
            Ordered as A, T, G, C. E.g: {'A': 0.22, 'T': 0.25, \
                'G': 0.25, 'C': 0.28}.
        mu_length (int): mean length of desired fragments
        std (int): standard deviation of desired fragment lengths

    Returns:
        list: list of selected terminal fragments
    """
    # calculated using https://www.nature.com/articles/srep04532#MOESM1
    term_frags = []
    for seq_id, seq in fasta.items():
        counts = seq_counts[seq_counts["seqID"] == seq_id]["count"]
        for _ in range(counts.iloc[0]):
            # pick no. of cuts from gauss fragment length distribution
            # non-uniformly random DNA fragmentation implementation based on
            # https://www.nature.com/articles/srep04532#Sec1
            # assume fragmentation by sonication for NGS workflow
            cuts = []
            cut_nucs = np.random.choice(list(nuc_probs.keys()),
                                        size=get_cut_number(len(seq),
                                                            mu_length),
                                        p=list(nuc_probs.values()))
            for nuc in cut_nucs:
                nuc_pos = [x.start() for x in re.finditer(nuc, seq)]
                pos = np.random.choice(nuc_pos)
                while pos in cuts:
                    pos = np.random.choice(nuc_pos)
                cuts.append(pos)

            cuts.sort()
            cuts.insert(0, 0)
            term_frag = ""

            # check if 3' fragment is in the correct size range
            if mu_length-std <= len(seq[cuts[-1]:len(seq)]) <= mu_length+std:
                term_frag = seq[cuts[-1]:len(seq)]

            if term_frag != "":
                term_frags.append(term_frag)

    return term_frags
