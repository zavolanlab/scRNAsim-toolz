
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
        for nt in range(len(seq),read_length):
            sequenced += choice(bases)
    else:
        for nt in range(read_length):
            sequenced += seq[nt]

    return sequenced

def simulate_sequencing(sequences, read_length):
    results = {}
    for index, key in enumerate(sequences):
        results[key] = read_sequence(sequences[key],read_length=read_length)

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

def run_read_sequencer(input_file_path, read_length, output_file_path):
    sequences = read_in_fasta(input_file_path)
    reads = simulate_sequencing(sequences, read_length)
    write_fasta(reads, output_file_path)
