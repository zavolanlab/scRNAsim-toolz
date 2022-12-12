import pandas as pd
from gtfparse import read_gtf

"""This script defines a BED from exon annotation in a GTF, to get sequences with transcript ID as header after usage in bedtools.

    For each transcript, take exons only and sort exons by start position (reverse order for -ve strand)
    Input: GTF file 
    Columns needed for BED: chr, start, end, transcript_id, score, strand, gene_id
    ...
    :returns: BED file format
    :rtype: dataframe
    """


gtf = read_gtf('../scrna-seq-simulation-main/inputs/ref_annotation.gtf')

gtf_exons = gtf[gtf["feature"] == "exon"]

gtf_exons = gtf_exons[["seqname", "start", "end", "transcript_id", "score", "strand", "gene_id"]]

gtf_df_neg = gtf_exons[gtf_exons["strand"] == "-"]
gtf_df_neg = gtf_df_neg.sort_values(['transcript_id','start'],ascending=False).groupby('transcript_id').head(len(gtf_df_neg. transcript_id))

gtf_df_pos = gtf_exons[gtf_exons["strand"] == "+"]
gtf_df_pos = gtf_df_pos.sort_values(['transcript_id','start'],ascending=True).groupby('transcript_id').head(len(gtf_df_pos. transcript_id))

pd.concat([gtf_df_pos, gtf_df_neg]).to_csv("bed_file.bed",sep="\t",index=False) #gtf_df_pos and gtf_df_neg must be dataframes

