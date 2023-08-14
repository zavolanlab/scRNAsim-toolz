"""
Changed on Tue Dec 20 14:06:20 2022

@author: RobinC
"""
import sys
from createprimer import CreatePrimer
from postprocessing import PostProcessRIBlast


def main():
    generate_RIBlast_input()

        

def generate_RIBlast_input():
    """This function creates a list of the filenames for the RIBlast"""
    my_primer = CreatePrimer()
    my_primer.create_fasta()
    primer_filename = my_primer.name +".fasta"
    transcripts_filename = "transcripts.fasta"
    
    return [primer_filename, transcripts_filename]

def create_gtf():
    gtf_file = PostProcessRIBlast().output
    print(gtf_file)

if __name__ == '__main__':
    main()
