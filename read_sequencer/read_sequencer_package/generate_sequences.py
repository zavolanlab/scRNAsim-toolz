import random

def generate_sequences(n, mean, sd):
    """
    Generates random sequences.

    Args:
        n (int): Amount of sequences to generate.
        mean (int): mean length of sequence (gaussian distribution).
        sd (float): standart deviation of length of sequence (gaussian distribution).

    Returns:
        list: of n sequences
    """
    dict1 = {}
    for i in range(n):
        keys = range(n)
        seq = ""
        nt = ["A", "T", "C", "G"]
        for value in range(round(random.gauss(mean, sd))):
            seq = seq + random.choice(nt)
        dict1[keys[i]] = seq
    return dict1

