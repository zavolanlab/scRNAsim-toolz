import argparse
import time
import transcript_sampler as ts

# exemple execution : python C:\...\final_exe.py  --input_gtf  "C:\...\input_files\test.gtf" --input_csv "C:\...\input_files\expression.csv"  --output_gtf "C:\...\output\output_gtf.gtf"  --output_csv "C:\...\ouput\output_gtf.gtf" --n_to_sample 100


def exe(input_gtf, input_csv, output_gtf, output_csv, transcript_nr, input_free=True):
    start = time.time()
    dict_repr_trans = ts.get_rep_trans(input_gtf)
    df_repr = ts.match_reprTranscript_expressionLevel(
        dict_reprTrans=dict_repr_trans, exprTrans=input_csv, gtf_file=input_gtf
    )
    print("Finiding match between representative transcripts and expression level file")
    print("Poisson sampling of transcripts")
    ts.transcript_sampling(transcript_nr, df_repr, output_csv)
    print("output csv file ready")
    print("writing output gtf file")
    ts.gtf_file_writer(input_gtf, dict_repr_trans, output_gtf)
    end = time.time()
    print("\nScript executed in {} sec\n".format(end - start))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="transcript sampler",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input_gtf", required=True, help="gtf file with genome annotation"
    )
    parser.add_argument(
        "--input_csv",
        required=True,
        help="csv or tsv file with transcript and their expression level ",
    )
    parser.add_argument(
        "--output_gtf",
        required=True,
        help="output path for the new gtf file of representative transcripts",
    )
    parser.add_argument(
        "--output_csv",
        required=True,
        help="output path for the new csv file of representative transcript and their sampled number",
    )
    parser.add_argument(
        "--n_to_sample", required=True, help="total number of transcripts to sample"
    )
    args = parser.parse_args()
    exe(
        args.input_gtf,
        args.input_csv,
        args.output_gtf,
        args.output_csv,
        args.n_to_sample,
    )
