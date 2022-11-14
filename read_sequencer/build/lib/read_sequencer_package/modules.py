def generate_sequences(n, mean, sd):
    """
    Generates random sequences.

    Args:
        n (int): Amount of sequences to generate.
        mean (int): mean length of sequence (gaussian distribution).
        sd (float): standard deviation of length of sequence (gaussian distribution).

    Returns:
        list: of n sequences
    """
    from random import gauss, choice
    dict = {}
    for i in range(n):
        keys = range(n)
        seq = ""
        nt = ["A", "T", "C", "G"]
        for value in range(abs(round(gauss(mean, sd)))):
            seq = seq + choice(nt)
        dict[keys[i]] = seq
    return dict


def read_in_fasta(file_path):
    '''
    This function reads in FASTA files.

    Args:
        file_path (str): A file path directing to the fasta file.  

    Returns:
        Dict: It returns a dictionary with sequences.

    '''
    sequences = {}
    f = open(file_path)
    for line in f:
        if line[0] == '>':
            defline = line.strip()
            defline = defline.replace('>', '')
        else:
            if defline not in sequences:
                sequences[defline] = ''
                sequences[defline] += line.strip()
    f.close()
    return sequences

def read_sequence(seq, read_length):
    '''
    This function reads a sequence of a specific read length and adds additional nucleotides if the sequence is 
    smaller then the requested length or cuts the sequence if its longer.

    Args:
        seq (str): the sequence to read 
        read_length (int): length of reads

    Returns:
        str: returns sequenced element

    '''
    from random import choice
    bases = ["A", "T", "C", "G"]
    sequenced = ''
    if read_length >= len(seq):
        for nt in range(len(seq)):
            sequenced += seq[nt]
        for nt in range(len(seq), read_length):
            sequenced += choice(bases)
    else:
        for nt in range(read_length):
            sequenced += seq[nt]

    return sequenced

def simulate_sequencing(sequences, read_length):
    """
    Simulates sequencing.

    Args:
        sequences (dict): Dictionary of sequences to sequence.
        read_length (int): length of reads

    Returns:
        dict: of n sequences as values 
    """
    results = {}
    for index, key in enumerate(sequences):
        results[key] = read_sequence(sequences[key], read_length=read_length)

    return results

import random
def generate_sequences(n, mean, sd):
    """
    Generates random sequences.

    Args:
        n (int): Amount of sequences to generate.
        mean (int): mean length of sequence (gaussian distribution).
        sd (float): standart deviation of length of sequence (gaussian distribution).

    Returns:
        dict: of n sequences
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

def write_fasta(sequences, file_path):
    """
    Takes a dictionary and writes it to a fasta file.
    Must specify the filename when calling the function.

    Args:
        sequences (dict): Dictionary of sequence.
        file_path (str): A file path directing to the output folder.
        
    """
    from textwrap import wrap
    with open(file_path, "w") as outfile:
        for key, value in sequences.items():
            outfile.write(key + "\n")
            outfile.write("\n".join(wrap(value, 60)))
            outfile.write("\n")

class read_sequencer:
    def __init__(self):
        self.sequences = {}
        self.reads = {}

    def add_random_sequences(self, n, mean, sd):
        self.sequences = generate_sequences(n, mean, sd)

    def read_fasta(self, input_file):
        self.sequences = read_in_fasta(input_file)

    def run_sequencing(self, read_length):
        self.reads = simulate_sequencing(self.sequences, read_length)

    def write_fasta(self, output_file_path):
        write_fasta(self.reads, output_file_path)
