# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:17:06 2022

@author: baerma
"""
from createprimer import CreatePrimer

def generate_RIBlast_input():
    """This function creates a list of the filenames for the RIBlast"""
    my_primer = CreatePrimer()
    my_primer.create_fasta()
    primer_filename = my_primer.name +".fasta"
    transcripts_filename = "transcripts.fasta"
    
    return [primer_filename, transcripts_filename]

print(generate_RIBlast_input())



