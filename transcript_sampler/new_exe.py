"""This module executes the transcript_sampler"""
import argparse
import time
import logging
logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO,
    )
from find_reptrans import FindRepTrans  # pylint: disable=E0401,C0413
from match_reptrans_explvl import MatchReptransExplvl  # pylint: disable=E0401,C0413
from poisson_sampling import SampleTranscript  # pylint: disable=E0401,C0413

find_rep_trans = FindRepTrans()
match_reptrs_explvl = MatchReptransExplvl()
poisson_sample = SampleTranscript()

LOG = logging.getLogger(__name__)


def exe(input_gtf, input_csv, output_gtf, output_csv, transcript_nr):
    """Execute transcript sampler."""
    start = time.time()
    LOG.info("Started transcript sampler...")
    dict_repr_trans = find_rep_trans.get_rep_trans(input_gtf)
    df_repr = match_reptrs_explvl.match_repr_transcript_expression_level(
        dict_reprTrans=dict_repr_trans, exprTrans=input_csv, gtf_file=input_gtf
        )
    LOG.info(
        "Finding match between representative transcripts \
            and expression level file"
        )
    LOG.info("Poisson sampling of transcripts")
    poisson_sample.transcript_sampling(transcript_nr, df_repr, output_csv)
    LOG.info("Output CSV file ready")

    LOG.info("Writing output GTF file")
    find_rep_trans.gtf_file_writer(input_gtf, dict_repr_trans, output_gtf)

    end = time.time()
    LOG.info("Script executed in %s sec", (end - start))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Transcript sampler",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input_gtf", required=True,
        help="GTF file with genome annotation"
        )
    parser.add_argument(
        "--input_csv", required=True,
        help="CSV or TSV file with transcripts and their expression level"
        )
    parser.add_argument(
        "--output_gtf", required=True,
        help="Output path for the new GTF file of representative transcripts"
        )
    parser.add_argument(
        "--output_csv", required=True,
        help="Output path for the new CSV file of representative transcripts \
            and their sampled number"
        )
    parser.add_argument(
        "--n_to_sample", required=True,
        help="Total number of transcripts to sample"
        )
    args = parser.parse_args()
    print(args)

    exe(
        args.input_gtf,
        args.input_csv,
        args.output_gtf,
        args.output_csv,
        args.n_to_sample,
    )
