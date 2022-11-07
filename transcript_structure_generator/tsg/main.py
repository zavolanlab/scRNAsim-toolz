import logging

import pandas as pd

LOG = logging.getLogger(__name__)


def read_abundances(transcripts_file: str) -> pd.DataFrame:
    """Read abundance file into dataframe

    Args:
        transcripts_file (str): Input filename

    Returns:
        pd.DataFrame: Transcript abundances ("id", "count")
    """
    cols: list = ["id", "count"]
    if transcripts_file.endswith(".tsv"):
        return pd.read_table(transcripts_file, header=None, names=cols)
    elif transcripts_file.endswith(".csv"):
        return pd.read_csv(transcripts_file, header=None, names=cols)


def filter_df(df: pd.DataFrame, transcripts: list = []) -> pd.DataFrame:
    # Filter annotations to exon and highest transcript support level.
    # If list of transcript ids is given, filter for that as well.
    df_filter = df[
        (df["feature"] == "exon")
        & (df["free_text"].str.contains('transcript_support_level "1"'))
    ]
    if len(transcripts) > 0:
        df_filter = df_filter.str.contains("|".join(transcripts), regex=True)

    return df_filter


def str_to_dict(s: str) -> dict:
    # split between key/value pairs
    # remove empty list items and split key, value pairs
    item_list: list = [x.split() for x in s.split(";") if len(x) > 0]
    # remove quotes for values and return dictionary
    return {item[0]: item[1].strip('"') for item in item_list}


def dict_to_str(d: dict) -> str:
    # join key, value pairs from dictionary with a space in a list,
    # then join items in list by ;
    # end on ;
    # check if value is nan
    s: str = (
        "; ".join([f'{key} "{value}"' for key, value in d.items() if value == value])
        + ";"
    )
    return s


def reverse_parse_free_text(df_all: pd.DataFrame) -> pd.DataFrame:
    # the first 8 columns should be constant according to gtf file standard
    # we assume that further columns are parsed free text columns
    df_free_text = df_all.iloc[:, 8:]

    df = df_all.iloc[:, :8]
    df["free_text"] = df_free_text.agg(pd.Series.to_dict, axis=1).apply(dict_to_str)
    return df


def write_gtf(df: pd.DataFrame, filename: str) -> None:
    # Make sure the data types are correct.
    df = df.astype(Gtf.dtypes)

    df.to_csv(
        filename,
        sep="\t",
        header=False,
        index=False,
        quoting=None,
        quotechar="'",
        mode="a",
    )


def write_header(annotations_file: str) -> None:
    with open(annotations_file, "w") as fh:
        fh.write("\t".join(Gtf.dtypes.keys()) + "\n")


class Gtf:
    """Class to read transcripts annotations file and parse it into a pandas Dataframe.

    Args:
        annotations_file: Path to gtf file.

    Attributes:
        annotations_file: File with transcript annotation of the genome

    """

    dtypes = {
        "seqname": object,
        "source": object,
        "feature": object,
        "start": int,
        "end": int,
        "score": object,
        "strand": object,
        "frame": object,
        "free_text": object,
    }

    def __init__(self):
        self.parsed = False
        self.original_columns = list(self.dtypes.keys())
        self.free_text_columns = []

    def read_file(self, annotations_file: str) -> None:
        # for large annotation files, iterate over lines and filter before saving to dataframe
        reader = pd.read_table(
            annotations_file,
            sep="\t",
            comment="#",
            names=self.dtypes.keys(),
            dtype=self.dtypes,
            chunksize=100000,
            iterator=True,
        )
        self.df = pd.concat([filter_df(chunk) for chunk in reader])

    def from_dataframe(df: pd.DataFrame) -> None:
        self.free_text_columns = [
            col for col in df.columns if col not in self.original_columns
        ]
        self.df = df
        if not "free_text" in df.columns:
            self.parsed = True

    def parse_free_text(self):
        assert self.parsed == False
        # create dataframe with columns for values in free_text column
        df_free_text = self.df["free_text"].map(str_to_dict).apply(pd.Series)
        # remember which columns come from free_text
        self.free_text_columns = df_free_text.columns
        # join free_text columns to original dataframe and drop the "free_text" column itself
        self.df = self.df.drop("free_text", axis=1)
        self.original_columns = self.df.columns
        self.df = self.df.join(df_free_text, how="inner")
        # remember that current dataframe is parsed, i.e. can't be written in gtf format
        self.parsed = True

    def reverse_parse_free_text(self):
        assert self.parsed == True
        # create dataframe with only free_text columns
        df_free_text = self.df[self.free_text_columns]
        # filter current dataframe to only original columns, except "free_text" column
        self.df = self.df[self.original_columns]
        # undo parsing and save result in "free_text" column
        self.df["free_text"] = df_free_text.agg(pd.Series.to_dict, axis=1).apply(
            dict_to_str
        )
        # remember that current dataframe is not parsed
        self.parsed = False

    def pick_transcript(self, transcript_id: str) -> pd.DataFrame:
        return self.df.query(f"transcript_id == '{transcript_id}'")

def sample_transcripts(
    input_transcripts_file: str,
    input_annotations_file: str,
    prob_inclusion: float,
    output_transcripts_file: str,
    output_annotations_file: str,
):
    transcripts = read_abundances(input_transcripts_file)

    annotations = Gtf()
    annotations.read_file(input_annotations_file)
    annotations.parse_free_text()
