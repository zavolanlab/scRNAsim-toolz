import argparse
import logging

from cdna import CDNAGen


def parser():
    parser = argparse.ArgumentParser(
        prog="cDNA generator",
        description="Generate cDNA sequences based on primer probabilities.",
    )
    parser.add_argument("-ifa", "--input_fasta", help="genome fasta file", required=True)
    parser.add_argument("-igtf", "--input_gtf", help="gtf file", required=True)
    parser.add_argument("-ofa", "--output_fasta", help="output fasta file", required=True)
    parser.add_argument("-icpn", "--input_copy_number", help="input copy number (csv) file", required=True)
    parser.add_argument("-ocsv", "--output_csv", help="output fasta file", required=True)
    args = parser.parse_args()

    print(' \n'.join(f'{k}={v}' for k, v in vars(args).items()))
    print()
    CDNA = CDNAGen(
        ifasta=args.input_fasta,
        igtf=args.input_gtf,
        icpn=args.input_copy_number,
        ocsv=args.output_csv,
        ofasta=args.output_fasta,
    )
    return CDNA


if __name__ == "__main__":
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO,
    )
    LOG = logging.getLogger(__name__)
    print("**********************")
    print("Running cDNA generator")
    print("**********************")
    cnda_object = parser()
