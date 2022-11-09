'''
This function reads in FASTA files

argument is file_path

it returns a dictionary with the sequences


'''


import sys

def read_in_fasta(file_path):
    sequences = {}
    f = open(file_path)
    for line in f:
        if line[0] == '>':
            defline = line.srtip()
            defline = defline.replace('>', '')
        else:
            if defline not in sequences:
                sequences[defline] = ''
            sequences[defline] += line.strip()
    return sequences       

