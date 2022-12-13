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
    LOG.info("sequence_extractor begins")
    fasta_list = exon_concatenation(args.input_fasta_file)
    final_list = poly_a_addition_to_fasta_list(fasta_list)
    with open(args.output_file_name, 'w') as fasta_out:
        fasta_out.write('\n'.join('%s\n%s' % x for x in final_list))
    LOG.info("sequence_extractor ends")

if ___name__ == 'main':
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO,
    )
    LOG = logging.getLogger(__name__)
    main()
