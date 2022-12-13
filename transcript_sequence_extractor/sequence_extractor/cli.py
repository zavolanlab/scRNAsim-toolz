""" command line script to be run on output fasta file from bedtools getfasta """
import argparse
import logging
from exon_concatenation import exon_concatenation
from poly_a import poly_a_addition_to_fasta_list

parser = argparse.ArgumentParser(
    prog = 'transcript_sequence_extractor',
    description = 'extracts transcript sequences from genome sequence and ouputs transcripts with PolyA tail added to them')
parser.add_argument('--input_fasta_file',
                    help='fasta file obtained from bedtools')
parser.add_argument('--output_file_name',
                    help='Name of the output fasta file')

args = parser.parse_args()

def main():
    """Runs on the output from bedtools and concatenates the exons together and adds a polyA tail and outputs a fasta file.

    Args:
        None: this will run on its own by taking the information from argparse

    Returns:
        A fasta file with a single entry for each transcript ID with polyA tail being added onto the sequence at 3'end
    """
    LOG.info("sequence_extractor begins")
    fasta_list = exon_concatenation(args.input_fasta_file)
    final_list = poly_a_addition_to_fasta_list(fasta_list)
    with open(args.output_file_name, 'w', encoding="utf-8") as fasta_out:
        fasta_out.write('\n'.join('%s\n%s' % x for x in final_list))
    LOG.info("sequence_extractor ends")

if ___name__ == 'main':
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO,
    )
    LOG = logging.getLogger(__name__)
    main()
