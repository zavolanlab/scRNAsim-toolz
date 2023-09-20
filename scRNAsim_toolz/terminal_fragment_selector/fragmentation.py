"""Fragment sequences."""
import random
import pandas as pd  # type: ignore


def fragmentation(fasta: dict, seq_counts: pd.DataFrame,
                  mean: int = 100, std: int = 10) -> list:
    """Fragment cDNA sequences and select terminal fragment.

    Args:
        fasta_file (dict): dictionary of {transcript IDs: sequences}
        counts_file (pd.DataFrame): dataframe with sequence counts and IDs
        mean (int): mean length of desired fragments
        std (int): standard deviation of desired fragment lengths

    Returns:
        list: list of selected terminal fragments
    """
    if not isinstance(mean, int):
        raise ValueError(f"Mean must be numeric, not {type(mean)}")

    if not isinstance(std, int):
        raise ValueError(f"Std must be numeric, not {type(mean)}")

    term_frags = []
    for seq_id, seq in fasta.items():
        counts = seq_counts[seq_counts["seqID"] == seq_id]["count"]
        for _ in range(counts.iloc[0]):
            cuts = []
            seq_len = len(seq)
            prob_cut_per_base = 1/mean

            for i in range(seq_len):
                rand_prob = random.uniform(0, 1)
                if rand_prob < prob_cut_per_base:
                    cuts.append(i)

            cuts.sort()
            cuts.insert(0, 0)
            term_frag = ""

            # check if 3' fragment is in the correct size range
            if mean - std <= len(seq[cuts[-1]:len(seq)]) <= mean + std:
                term_frag = seq[cuts[-1]:len(seq)]

            if term_frag != "":
                term_frags.append(term_frag)

    return term_frags
