"""Sample transcripts by Poisson-sampling."""

import pandas as pd  # type: ignore
import numpy as np  # type: ignore


# pylint: disable=R0903
class SampleTranscript:
    """Sample transcript.

    This part of the code does Poisson sampling proportionally
    to gene expression levels for each gene.

    input:  total transcript number (int)
            csv file with gene id and  gene expression levels
            (columns named 'id' and 'level')

    output: csv file with gene id and count
            gtf file with transcript samples
    """

    @staticmethod
    def transcript_sampling(total_transcript_number, df_repr, output_csv):
        """Sample transcript based on Poisson-sampling."""
        total = df_repr["level"].sum()
        total_transcript_number = int(total_transcript_number)
        normalized = total_transcript_number / total
        levels = np.random.poisson(df_repr["level"] * normalized)
        transcript_numbers = pd.DataFrame({
            "id": df_repr["id"], "count": levels
            })
        transcript_numbers.to_csv(output_csv, index=False, header=False)
