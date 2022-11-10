
import random

def generate_sequences(n, mean, sd):
    """Summary line.

    Generates random sequences.

    Args:
        n (int): Amount of sequences to generate.
        mean (int): mean length of sequence (gaussian distribution).
        sd (float): standart deviation of length of sequence (gaussian distribution).

    Returns:
        list: of n sequences
    """
    l1 = []
    for i in range(n):
        seq = ""
        nt = ["A", "T", "C", "G"]
        for pos in range(round(random.gauss(mean, sd))):
            seq = seq + random.choice(nt)
        l1.append(seq)
    return l1

