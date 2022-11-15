import argparse
import logging
from pre_bedtools import exon_extraction_from_gtf
from exon_concatenation import exon_concatenation
from polyA import PolyA_generator
from list_to_file import list_to_file

parser = argparse.ArgumentParser(
    prog = 'transcript_sequence_extractor',
    description = 'extracts transcript sequences from genome sequence and ouputs transcripts with PolyA tail added to them')
parser.add_argument('--input_fasta_file',
                    help='genome fasta file')
parser.add_argument('--input_gtf',
                    help='gtf file')
parser.add_argument('--output_file_name',
                    help='output fasta file')

args = parser.parse_args()

def main():
    LOG.info("sequence_extractor begins")
    exon_extraction_from_gtf()
    exon_concatenation()
    PolyA_generator()
    list_to_file()
    LOG.info("sequence_extractor ends")

if ___name__ == 'main':
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO,
    )
    LOG = logging.getLogger(__name__)
    main()
