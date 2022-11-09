'''
This function reads sequences
arguments: seq is a list of sequences 
padding_probabilities is a number??

returns sequenced element

'''

import random



def read_sequence(seq, padding_probabilities, read_length):
    reading_element = random.choice(seq)

    bases =["A", "T", "C", "G"]

    if read_length > len(reading_element):
        for nt in [0:len(reading_element)]:
            sequenced += reading_element[nt]

        for nt2 in [len(reading_element):read_length]:
            sequenced += random.choice(bases)

    else: 
        for nt in [0:read_length]
        sequenced += reading_element[nt]


    return sequenced