"""Tests for main module."""

import pandas as pd
from tsg.main import Gtf, TranscriptGenerator, dict_to_str, str_to_dict


class TestFreeTextParsing:
    """Test if free text dictionary is correctly parsed."""

    def test_str2dict(self):
        """Test for str2dict function."""
        res = str_to_dict(
            'gene_id "GENE2"; transcript_id "TRANSCRIPT2"; \
                exon_number "1"; exon_id "EXON1";'
        )

        assert res == {
            "gene_id": "GENE2",
            "transcript_id": "TRANSCRIPT2",
            "exon_number": "1",
            "exon_id": "EXON1",
        }

    def test_dict2str(self):
        """Test for dict2str function."""
        res = dict_to_str(
            {
                "gene_id": "GENE2",
                "transcript_id": "TRANSCRIPT2",
                "exon_number": "1",
                "exon_id": "EXON1",
            }
        )
        print(res)
        assert res == (
            'gene_id "GENE2"; '
            'transcript_id "TRANSCRIPT2"; '
            'exon_number "1"; '
            'exon_id "EXON1";'
        )


class TestGtf:
    """Test if Gtf class works correctly."""

    cols = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "free_text",
    ]

    def test_init(self):
        """Test for init function."""
        annotations = Gtf()
        annotations.read_file("tests/resources/Annotation1.gtf")

        assert annotations.parsed is False
        assert annotations.original_columns == self.cols
        assert annotations.free_text_columns == []

    def test_parsed(self):
        """Test for parsed function."""
        annotations = Gtf()
        annotations.read_file("tests/resources/Annotation1.gtf")
        annotations.parse_key_value()

        assert annotations.parsed is True
        assert set(annotations.free_text_columns) == set(
            [
                "gene_id",
                "transcript_id",
                "exon_number",
                "exon_id",
                "transcript_support_level",
            ]
        )
        assert set(annotations.original_columns) == set(
            ["seqname", "source", "feature", "start",
             "end", "score", "strand", "frame"]
        )


class TestTranscriptGenerator:
    """Test for TranscriptGenerator class."""

    cols = [
        "start",
        "end",
        "strand",
        "transcript_id",
    ]

    df1 = pd.DataFrame(
        {
            "start": [1, 50, 80],
            "end": [20, 70, 100],
            "strand": ["+", "+", "+"],
            "exon_id": ["EXON1", "EXON2", "EXON3"],
        }
    )
    df2 = pd.DataFrame(columns=["start", "end", "strand"])

    def test_init(self):
        """Test for init."""
        transcripts = TranscriptGenerator("TRANSCRIPT1", 3, self.df1, 0.05)

        assert transcripts.strand == "+"

    def test_inclusions(self):
        """Test for inclusions."""
        transcripts = TranscriptGenerator("TRANSCRIPT1", 3, self.df1, 0.5)
        res = transcripts.get_inclusions()

        assert res.shape == (3, 3)

    def test_unique_inclusions(self):
        """Test for unique inclusions."""
        transcripts = TranscriptGenerator("TRANSCRIPT1", 3, self.df1, 0.5)
        transcripts.get_unique_inclusions()

    def test_get_df(self):
        """Test for get_df function."""
        inclusions = [False, True, False]
        expected_end = pd.Series([20, 79, 100], name="end")
        transcript_id = "TRANSCRIPT1_1"

        transcripts = TranscriptGenerator("TRANSCRIPT1", 3, self.df1, 0.5)
        res = transcripts.get_df(inclusions, transcript_id)

        assert res["transcript_id"].unique().item() == "TRANSCRIPT1_1"
        assert res["strand"].unique().item() == "+"
        assert res["exon_id"].tolist() == ["EXON1", "EXON2_1", "EXON3"]
        pd.testing.assert_series_equal(res["start"], self.df1["start"])
        pd.testing.assert_series_equal(res["end"], expected_end)
