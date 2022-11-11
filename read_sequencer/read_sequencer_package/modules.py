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
    This function reads in FASTA files

    argument is file_path

    it returns a dictionary with the sequences

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

def read_sequence(seq, read_length, padding_probabilities=None):
    '''
    This function reads sequences
    arguments: seq is a list of sequences
    padding_probabilities is a number??

    returns sequenced element

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
    results = {}
    for index, key in enumerate(sequences):
        results[key] = read_sequence(sequences[key], read_length=read_length)

    return results

def write_fasta(sequences, file_path):
    """
    Takes a dictionary and writes it to a fasta file
    Must specify the filename when caling the function
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
