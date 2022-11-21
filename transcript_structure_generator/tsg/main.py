"""Sample transcripts."""

import logging

import numpy as np
import pandas as pd
from tqdm import tqdm


LOG = logging.getLogger(__name__)


def read_abundances(transcripts_file: str) -> pd.DataFrame:
    """Read transcript-abundance file into dataframe.

    Args:
        transcripts_file (str): Input filename

    Returns:
        pd.DataFrame: Transcript abundances ("id", "count")

    Raises:
        ValueError: When the input file is neither csv or tsv
    """
    cols: list = ["id", "count"]
    if transcripts_file.endswith(".tsv"):
        return pd.read_table(transcripts_file, header=None, names=cols)
    elif transcripts_file.endswith(".csv"):
        return pd.read_csv(transcripts_file, header=None, names=cols)
    else:
        raise ValueError("File type needs to be either csv or tsv")


def filter_df(df: pd.DataFrame, transcripts: list = []) -> pd.DataFrame:
    """Filter annotations to include only exons with the highest transcript support level, i.e. TSL1.

    `feature` column is filtered on value "exon" and
    `free_text` column is filtered to include the string denoting the highest transcript support level
    ('transcript_support_level "1"').

    If a list of transcript IDs is given, `free_text` column is filtered to include one of the IDs.

    Args:
        df: A pd.DataFrame containing an unparsed gtf-file
        transcript: list of transcript IDs

    Returns:
        A pd.DataFrame containing only rows with exon annotations of highest transcript support level and,
        if provided, belonging to one of the given transcripts
    """
    df_filter = df[
        (df["feature"] == "exon")
        & (df["free_text"].str.contains('transcript_support_level "1"'))
    ]
    if len(transcripts) > 0:
        df_filter = df_filter["free_text"].str.contains(
            "|".join(transcripts), regex=True
        )

    return df_filter


def str_to_dict(s: str) -> dict:
    """Split between key/value pairs.

    Split string based on delimiter ';' into items, remove empty items and split items on delimiter ' ' into key/value pairs.
    Remove quotes from value strings and create a dictionary.

    Args:
        s: A string of the form 'gene_id "GENE1"; transcript_id "TRANSCRIPT1";'

    Returns:
        A dictionary containing e.g. {'gene_id': 'GENE1', 'transcript_id': 'TRANSCRIPT1'}
    """
    # split into items
    # remove empty items
    # split items into key/value pairs
    item_list: list = [x.split() for x in s.split(";") if len(x) > 0]
    # remove quotes for values and return dictionary
    return {item[0]: item[1].strip('"') for item in item_list}


def dict_to_str(d: dict) -> str:
    """Parse dictionary in gtf free_text column format.

    Takes e.g. dictionary {'gene_id': 'GENE1', 'transcript_id': 'TRANSCRIPT1'} and returns
    string 'gene_id "GENE1"; transcript_id "TRANSCRIPT1";'.
    Key/value pairs are joined by space to form an item and items are joinded by ';' to form a string.
    If a value is Not a Number (nan), the key/value pair is omitted from the string.

    Args:
        d: A dictionary of the form {'gene_id': 'GENE1', 'transcript_id': 'TRANSCRIPT1'}

    Returns:
        A string, e.g. 'gene_id "GENE1"; transcript_id "TRANSCRIPT1";'.
    """
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
    """Reverse parsing of gtf based pd.DataFrame to include only columns that are well defnined by gtf file standards.

    The first 8 defined columns are constant as defined by gtf file standards.
    Further columns are assumed to be parsed free-text columns (see Gtf.parse_free_text()).
    The parsed free-text columns are aggregated as a dictionary and the dictionry is parsed as a string in gtf format.

    Args:
        df_all: A pd.DataFrame containing a parsed gtf file.

    Returns:
        A DataFrame with the columns as defined by gtf file standards.
    """
    # Define pd.DataFrame containing only parsed free-text columns
    df_free_text = df_all.iloc[:, 8:]
    # Define pd.DataFrame containing only non-parsed columns
    df = df_all.iloc[:, :8]
    # Reverse parsing of free-text columns and add the result as column `free_text` to output pd.DataFrame
    df["free_text"] = df_free_text.agg(pd.Series.to_dict, axis=1).apply(dict_to_str)
    return df


def write_gtf(df: pd.DataFrame, filename: str) -> None:
    """Save a Gtf object to file in gtf format.

    Makes sure data types are correct and saves object in gtf format.

    Args:
        df: A pd.DataFrame containing a gtf file.
        Filename: File to save to.
    """
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
    """Write the header of an annotations file, consisting of the tab delimited column names.

    Args:
        annotations_file: Filename to write header to.
    """
    with open(annotations_file, "w") as fh:
        fh.write("\t".join(Gtf.dtypes.keys()) + "\n")


class Gtf:
    """Class to read transcripts annotations file into a Gtf object.

    Attributes:
        dtypes: A dictionary containing column names and respective data types.
        parsed: A boolean indicating if the pd.DataFrame is parsed.
        original_columns: A list of columns not touched by parsing.
        free_text_columns: A list of columns created during parsing of column `free_text`.
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
        """Initialize Gtf object."""
        self.parsed = False
        self.original_columns = list(self.dtypes.keys())
        self.free_text_columns = []

    def read_file(self, annotations_file: str) -> None:
        """Read gtf file.

        Iterate over chunks of the gtf file reading 100000 rows at a time. Filter chunks for exon annotations of
        the highest transcript support level. Concatenate chunks to get resulting pd.DataFrame.

        Args:
            annotations_file: Filename of annotations.

        Raises:
            ValueError: The file type is required to be gtf.
        """
        if not annotations_file.endswith("gtf"):
            raise ValueError("File type needs to be gtf")

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

    def from_dataframe(self, df: pd.DataFrame) -> None:
        """Initialize Gtf object from pandas Dataframe.

        Part of initialization is:
        Set dataframe attribute
        Check which columns belong to the free-text part of the GTF-file.
        Check if there are no columns called free-text and if so, sets the value of parsed attribute to TRUE.

        Args:
            df: pd.DataFrame
        """
        self.free_text_columns = [
            col for col in df.columns if col not in self.original_columns
        ]
        self.df = df
        if "free_text" not in df.columns:
            self.parsed = True

    def parse_free_text(self):
        """Parse key/value pairs from `free_text` column into column `key` with row entry `value`.

        Creates a dataframe with columns for keys in the free-text column instead of `free_text` column.
        Saves it to Gtf.df attribute.
        """
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
        """Reverses parsing of `free_text` column.

        Creates a data frame that can be written in gtf format to file. Parsed free-text columns are aggregated
        into `free_text` column according to gtf format specification.
        """
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
        """Filter annotations to a given transcript ID."""
        return self.df.query(f"transcript_id == '{transcript_id}'")


class TranscriptGenerator:
    """Class to sample a transcript."""

    def __init__(
        self,
        transcript_id: str,
        transcript_count: int,
        transcript_df: pd.DataFrame,
        prob_inclusion: float,
    ):
        """Initialize TranscriptGenerator object."""
        assert len(transcript_df) > 0
        assert transcript_count > 0
        assert (prob_inclusion >= 0) and (prob_inclusion <= 1)

        self.id = transcript_id
        self.count = transcript_count
        self.df = transcript_df
        self.no_exons = len(transcript_df)
        self.strand = self.df["strand"].unique().item()
        self.prob_inclusion = prob_inclusion

    def _get_inclusions(self) -> np.array:
        """Generate inclusions array.

        Each column corresponds to one sample and the number of columns corresponds to the number of samples.

        Returns:
            inclusion_arr: A boolean np.array, where True means intron inclusion.
        """
        inclusion_arr = np.random.rand(self.no_exons, self.count) < self.prob_inclusion
        if self.strand == "+":
            inclusion_arr[-1, :] = False
        elif self.strand == "-":
            inclusion_arr[-1, :] = False

        return inclusion_arr

    def _get_unique_inclusions(self) -> tuple[list, np.array, np.array]:
        """Inclusion of unique intron inclusion via arrays and counts and name generation of each unique count.

        Args:

        Returns:
            names: List of names for generated exons.
            inclusion_arr_unique: A boolean np.array where columns correspond to generated transcripts and rows to
                intron inclusion.
            counts: A np.array containing sample number per generated inclusions, i.e. transcript.
        """
        inclusion_arr = self._get_inclusions()
        # Unique intron inclusion arrays and counts
        inclusion_arr_unique, counts = np.unique(
            inclusion_arr, axis=1, return_counts=True
        )
        # Name for each generated transcript
        names = []
        for i in range(inclusion_arr_unique.shape[1]):
            if np.all(inclusion_arr_unique[:, i] == False, axis=0):
                names.append(self.id)
            else:
                names.append(f"{self.id}_{i}")

        return names, inclusion_arr_unique, counts

    def _get_df(self, inclusions: np.array, transcript_id: str) -> pd.DataFrame:
        """Take as input a dataframe filtered to one transcript and a boolean vector denoting intron inclusions.

        Args:
            inclusions (np.array): A boolean vector denoting intron inclusion.
            transcript_id (str): The transcript id.

        Returns:
            The generated transcript as a pd.DataFrame.
        """
        df_generated = self.df.copy()
        if self.strand == "+":
            origninal_end = df_generated["end"]
            df_generated["end"] = np.where(
                inclusions,
                df_generated["start"].shift(periods=-1, fill_value=-1) - 1,
                origninal_end,
            )
        if self.strand == "-":
            origninal_start = df_generated["start"]
            df_generated["start"] = np.where(
                inclusions,
                df_generated["end"].shift(periods=-1, fill_value=-1) + 1,
                origninal_start,
            )

        original_id = df_generated["exon_id"]
        df_generated["exon_id"] = np.where(
            inclusions,
            df_generated["exon_id"] + "_" + np.arange(len(df_generated)).astype(str),
            original_id,
        )

        df_generated["transcript_id"] = transcript_id
        return df_generated

    def generate_transcripts(self, filename: str) -> None:
        """Write transcripts to file.

        Args:
            filename (str): Output csv filename.
        """
        ids, inclusions, counts = self._get_unique_inclusions()
        with open(filename, "a") as fh:
            for transcript_id, transcript_count in zip(ids, counts):
                fh.write(f"{transcript_id},{self.id},{transcript_count}\n")

    def generate_annotations(self, filename: str) -> None:
        """Generate a annotations in gtf format for sampled transcript.

        Args:
            Filename: Output gtf filename.

        Raises:
            ValueError: If given transcript ID could not be sampled.
        """
        ids, inclusions, counts = self._get_unique_inclusions()
        n_unique = len(ids)

        try:
            df = pd.concat(
                [self._get_df(inclusions[:, i], ids[i]) for i in range(n_unique)]
            )
            df = reverse_parse_free_text(df)

            write_gtf(df, filename)
            LOG.debug(f"Transcript {self.id} sampled")
        except ValueError:
            LOG.error(f"Transcript {self.id} could not be sampled.")


def sample_transcripts(
    input_transcripts_file: str,
    input_annotations_file: str,
    prob_inclusion: float,
    output_transcripts_file: str,
    output_annotations_file: str,
):
    """Read input files, iterate over transcript IDs, sample each transcript and save results.

    Args:
        input_transcripts_file (str): Filename of transcript abundances, needs to be csv or tsv.
        input_annotations_file (str): Filename of annotations, needs to be gtf.
        prob_inclusion (float): Probability of intron inclusion, needs to be float in range [0,1].
        output_transcripts_file (str): Filename of file to write sampled transcripts to.
        output_annotations_file (str): Filename of file to write generated annotations to.
    """
    transcripts = read_abundances(input_transcripts_file)

    LOG.info("Parsing annotations...")
    annotations = Gtf()
    annotations.read_file(input_annotations_file)
    annotations.parse_free_text()
    LOG.info("Done parsing...")

    LOG.info("Start sampling transcripts...")
    # Set up output file, write header once and append data in loop
    write_header(output_annotations_file)

    for _, row in tqdm(transcripts.iterrows()):
        transcript_id = row["id"]
        transcript_count = row["count"]

        transcript_df = annotations.pick_transcript(transcript_id)
        transcripts = TranscriptGenerator(
            transcript_id,
            transcript_count,
            transcript_df,
            prob_inclusion=prob_inclusion,
        )
        transcripts.generate_annotations(output_annotations_file)
        transcripts.generate_transcripts(output_transcripts_file)
    LOG.info("Done.")
