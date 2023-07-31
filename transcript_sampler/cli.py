"""This module executes the transcript_sampler."""
import argparse
import time
import logging

logging.basicConfig(
    format='[%(asctime)s: %(levelname)s] %(message)s \
        (module "%(module)s")',
    level=logging.INFO,
    )

from transcript_sampler.find_reptrans import FindRepTrans  # noqa: E402,E501 # pylint:disable=wrong-import-position
from transcript_sampler.match_reptrans_explvl import MatchReptransExplvl  # noqa: E402,E501 # pylint:disable=wrong-import-position
from transcript_sampler.poisson_sampling import SampleTranscript  # noqa: E402,E501 # pylint:disable=wrong-import-position

find_rep_trans = FindRepTrans()
match_reptrs_explvl = MatchReptransExplvl()
poisson_sample = SampleTranscript()


def main():
    """Execute transcript sampler."""
    parser = argparse.ArgumentParser(
        description="Transcript sampler",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input_gtf", required=True, default=None,
        help="GTF file with genome annotation"
        )
    parser.add_argument(
        "--input_csv", required=True, default=None,
        help="CSV or TSV file with transcripts and their expression level"
        )
    parser.add_argument(
        "--output_gtf", required=True, default=None,
        help="Output path for the new GTF file of representative transcripts"
        )
    parser.add_argument(
        "--output_csv", required=True, default=None,
        help="Output path for the new CSV file of representative transcripts \
            and their sampled number"
        )
    parser.add_argument(
        "--n_to_sample", required=True, default=None,
        help="Total number of transcripts to sample"
        )
    args = parser.parse_args()

    log = logging.getLogger("main")
    start = time.time()
    log.info("Started transcript sampler...")
    dict_repr_trans = find_rep_trans.get_rep_trans(args.input_gtf)
    df_repr = match_reptrs_explvl.match_repr_transcript_expression_level(
        dict_reprTrans=dict_repr_trans,
        exprTrans=args.input_csv,
        gtf_file=args.input_gtf
        )
    log.info(
        "Finding match between representative transcripts \
            and expression level file"
        )
    log.info("Poisson sampling of transcripts")
    poisson_sample.transcript_sampling(
        args.n_to_sample, df_repr, args.output_csv)
    log.info("Output CSV file ready")

    log.info("Writing output GTF file")
    find_rep_trans.gtf_file_writer(
        args.input_gtf, dict_repr_trans, args.output_gtf)

    end = time.time()
    log.info("Script executed in %s sec", (end - start))


if __name__ == "__main__":
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s \
            (module "%(module)s")',
        level=logging.INFO,
        )
    main()
