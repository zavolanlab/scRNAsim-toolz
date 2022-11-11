import pandas as pd
import numpy as np
import argparse


'''
Sample transcript 

This part of the code does Poisson sampling proportionally to gene expression levels for each gene. 
 
input:  total transcript number (int) 
        csv file with gene id and  gene expression levels (columns named 'id' and 'level')

output: csv file with gene id and count
        gtf file with transcript samples
'''


def transcript_sampling(total_transcript_number, csv_file, output_csv):
    df = pd.read_csv(csv_file, sep='\t', lineterminator='\n', names=["id", "level"])
    levels = []
    sums = df['level'].tolist()
    total = sum(sums)
    normalized = total_transcript_number/total
    for expression_level in df['level']:
        poisson_sampled = np.random.poisson(expression_level*normalized)
        levels.append(poisson_sampled)

    transcript_numbers = pd.DataFrame({'id': df['id'],'count': levels})
    pd.DataFrame.to_csv(transcript_numbers, output_csv)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Transcript Poisson sampler, csv output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--expression_level", required=True, help="csv file with expression level")
    parser.add_argument("--output_csv", required=True, help="output csv file")
    parser.add_argument("--input_csv", required=True, help="input csv file")
    parser.add_argument("--transcript_number", required=True, help="total number of transcripts to sample")
    args = parser.parse_args()


    transcript_sampling(args.transcript_number, args.input_csv, args.output_csv, args.transcript_number)


