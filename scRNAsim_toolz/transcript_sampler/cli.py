"""This module executes the transcript_sampler."""
import argparse
import time
import logging
from scRNAsim_toolz.version import __version__

logging.basicConfig(
    format='[%(asctime)s: %(levelname)s] %(message)s \
        (module "%(module)s")',
    level=logging.INFO,
    )

from .find_reptrans import FindRepTrans  # noqa: E402,E501 # pylint:disable=wrong-import-position # type: ignore
from .match_explvl import MatchReptransExplvl  # noqa: E402,E501 # pylint:disable=wrong-import-position
from .poisson_sampling import SampleTranscript  # noqa: E402,E501 # pylint:disable=wrong-import-position

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
        "-ic", "--input_csv", required=True, default=None,
        help="CSV or TSV file with transcripts and their expression level"
        )
    parser.add_argument(
        "-ig", "--input_gtf", required=True, default=None,
        help="GTF file with genome annotation"
        )
    parser.add_argument(
        "-oc", "--output_csv", required=True, default=None,
        help="Output path for the new CSV file of representative transcripts "
             "and their sampled number"
        )
    parser.add_argument(
        "-og", "--output_gtf", required=True, default=None,
        help="Output path for the new GTF file of representative transcripts"
        )
    parser.add_argument(
        "-n", "--n_to_sample", required=True, default=None,
        help="Total number of transcripts to sample"
        )
    parser.add_argument(
        '-v', '--version', action='version',
        version=f'scRNAsim version: {__version__}'
    )
    args = parser.parse_args()

    log = logging.getLogger("main")
    start = time.time()
    log.info("Started transcript sampler.")
    dict_repr_trans = find_rep_trans.get_rep_trans(args.input_gtf)
    df_repr = match_reptrs_explvl.match_repr_transcript_expression_level(
        dict_repr_trans=dict_repr_trans,
        expr_trans=args.input_csv,
        gtf_file=args.input_gtf
        )
    log.info(
        "Finding match between representative transcripts "
        "and expression level file..."
        )
    log.info("Poisson sampling of transcripts...")
    poisson_sample.transcript_sampling(
        args.n_to_sample, df_repr, args.output_csv)
    log.info("Output CSV file ready.")

    log.info("Writing output GTF file...")
    find_rep_trans.gtf_file_writer(
        args.input_gtf, dict_repr_trans, args.output_gtf)

    end = time.time()
    log.info("Script executed in %s sec.", round(end - start, 2))


if __name__ == "__main__":
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s \
            (module "%(module)s")',
        level=logging.INFO,
        )
    main()
