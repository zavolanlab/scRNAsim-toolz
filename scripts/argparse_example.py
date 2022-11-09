import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="transcript sampler",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--annotation", required=True, help="gtf file with genome annotation")
    parser.add_argument("--expression_level", required=True, help="csv file with expression level")
    parser.add_argument("--output_csv", required=True, help="output csv file")
    parser.add_argument("--output_gtf", required=True, help="output gtf file")
    parser.add_argument("--transcript_number", required=True, help="total number of transcripts to sample")
    args = parser.parse_args()



# script name being whatever we call our workflow class/function that takes the inputs
    script_name(args.annotation, args.expression_level, args.output_csv, args.output_gtf, args.transcript_number)


"""note:
so to run our entire workflow, it will suffice to type the following into the command line:
        python filename.py \
            --annotation filename1.gtf \
            --expression_level filename2.gtf \
            --output_csv output1.csv \
            --output_gtf output2.gtf \
            --transcript_number 100 (some number)


"""