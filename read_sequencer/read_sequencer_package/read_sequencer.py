import logging
from random import choice, gauss
from textwrap import wrap

LOG = logging.getLogger(__name__)


def read_in_fasta(file_path: str) -> dict[str,str]:
    """
    This function reads in FASTA files.

    Args:
        file_path: A file path directing to the fasta file.

    Returns:
        Dict: It returns a dictionary with sequences.
    """
    LOG.info("Reading in FASTA files from destination.")
    sequences: dict[str,str] = {}
    f = open(file_path)
    for line in f:
        if line[0] == '>':
            def_line = line.strip()
            def_line = def_line.replace('>', '')
        else:
            if def_line not in sequences:
                sequences[def_line] = ''
                sequences[def_line] += line.strip()
    f.close()
    return sequences

def read_sequence(seq:str, read_length:int) -> str:
    """
    This function reads a sequence of a specific read length and adds additional nucleotides if the sequence is 
    smaller than the requested length or cuts the sequence if its longer.

    Args:
        seq: the sequence to read
        read_length: length of reads

    Returns:
        str: returns sequenced element

    """
    bases: list[str] = ["A", "T", "C", "G"]
    sequenced: str = ''
    if read_length > len(seq):
        for nt in range(len(seq)):
            sequenced += seq[nt]
        for nt in range(len(seq), read_length):
            sequenced += choice(bases)
    else:
        for nt in range(read_length):
            sequenced += seq[nt]

    return sequenced

def simulate_sequencing(sequences: dict[str,str], read_length: int) -> dict[str,str]:
    """
    Simulates sequencing.

    Args:
        sequences: Dictionary of sequences to sequence.
        read_length: length of reads

    Returns:
        dict: of n sequences as values 
    """
    LOG.info("Sequencing in progress....")
    results: dict[str,str] = {}
    for index, key in enumerate(sequences):
        results[key] = read_sequence(sequences[key], read_length=read_length)
    LOG.info("Sequencing was successfully executed.")
    return results

def generate_sequences(n: int, mean: int, sd: int) -> dict[str,str]:
    """
    Generates random sequences.

    Args:
        n: Amount of sequences to generate.
        mean: mean length of sequence (gaussian distribution).
        sd: standard deviation of length of sequence (gaussian distribution).

    Returns:
        dict: of n sequences
    """
    LOG.info("Generating random sequences.")
    sequences: dict[str,str] = {}
    for i in range(n):
        seq: str = ""
        bases: list[str] = ["A", "T", "C", "G"]
        for nt in range(abs(round(gauss(mean, sd)))):
            seq = seq + choice(bases)
        key: str = str(i) + ': length ' + str(len(seq)) + ' nt'
        sequences[key] = seq
    return sequences

def write_fasta(sequences: dict[str,str], file_path: str) -> None:
    """
    Takes a dictionary and writes it to a fasta file.
    Must specify the filename when calling the function.

    Args:
        sequences: Dictionary of sequence.
        file_path: A file path directing to the output folder.
        
    """
    LOG.info("Writing FASTA file.")
    with open(file_path, "w") as outfile:
        for key, value in sequences.items():
            outfile.write(key + "\n")
            outfile.write("\n".join(wrap(value, 60)))
            outfile.write("\n")

class ReadSequencer:
    def __init__(self):
        self.sequences: dict[str,str] = {}
        self.reads: dict[str,str] = {}

    def add_random_sequences(self, n: int, mean: int, sd: int):
        self.sequences: dict[str,str] = generate_sequences(n, mean, sd)

    def read_fasta(self, input_file):
        self.sequences: dict[str,str] = read_in_fasta(input_file)

    def run_sequencing(self, read_length: int):
        self.reads: dict[str,str] = simulate_sequencing(self.sequences, read_length)

    def write_fasta(self, output_file_path: str):
        write_fasta(self.reads, output_file_path)
