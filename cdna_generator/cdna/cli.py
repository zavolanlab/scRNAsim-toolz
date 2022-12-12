import cdna
import argparse


def parser():
    parser = argparse.ArgumentParser(
        prog="cDNA generator",
        description="Generate cDNA sequences based on primer probabilities.",
    )
    parser.add_argument("--input_fasta_file", help="genome fasta file")
    parser.add_argument("--input_gtf", help="gtf file")
    parser.add_argument("--output_fasta_name", help="output fasta file")
    parser.add_argument("--input_copy_number", help="input copy number (csv) file")
    parser.add_argument("--output_csv_name", help="output fasta file")
    args = parser.parse_args()
    CDNA = cdna.cdna.CDNAGen(
        ifasta=args["input_fasta_file"],
        igtf=args["input_gtf_file"],
        icpn=args["input_copy_number"],
        ocsv=args["output_csv_name"],
        ofasta=args["output_fasta_name"],
    )
    return CDNA
