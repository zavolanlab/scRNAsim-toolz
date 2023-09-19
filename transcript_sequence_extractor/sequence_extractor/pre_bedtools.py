"""Prepare gtf file for bedtools.

This script defines a BED from exon annotation in a GTF, to get exon
coordinates for use in bedtools. It also ensures that the concatenation happens
in the correct order, regardless of the strandedness of the transcript.

Args:
    GTF file

Returns:
    BED file with the format:
        chr, start, end, transcript_id, score, strand, gene_id
"""
import pandas as pd  # type: ignore
from gtfparse import read_gtf  # type: ignore


def pre_bedtools_mode(args):
    """Execute the 'pre_bedtools' mode."""
    gtf = read_gtf(args.input_gtf_file, result_type="pandas")
    gtf_exons = gtf[gtf["feature"] == "exon"]
    gtf_exons = gtf_exons[
        ["seqname", "start", "end",
         "transcript_id", "score",
         "strand", "gene_id"]
    ]

    gtf_df_neg = gtf_exons[gtf_exons["strand"] == "-"]
    gtf_df_neg = (
        gtf_df_neg.sort_values(["transcript_id", "start"], ascending=False)
        .groupby("transcript_id")
        .head(len(gtf_df_neg.transcript_id))
    )

    gtf_df_pos = gtf_exons[gtf_exons["strand"] == "+"]
    gtf_df_pos = (
        gtf_df_pos.sort_values(["transcript_id", "start"], ascending=True)
        .groupby("transcript_id")
        .head(len(gtf_df_pos.transcript_id))
    )

    if args.output_bed_file is None:
        args.output_bed_file = "default_output.bed"

    pd.concat([gtf_df_pos, gtf_df_neg]).to_csv(
        args.output_bed_file, sep="\t", index=False
    )  # gtf_df_pos and gtf_df_neg must be dataframes
