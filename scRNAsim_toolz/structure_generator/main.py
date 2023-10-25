"""Sample transcripts."""

import logging

import numpy as np  # type: ignore
import pandas as pd  # type: ignore

LOG = logging.getLogger(__name__)


def read_abundances(transcripts_file: str) -> pd.DataFrame:
    """Read transcript-abundance file into dataframe.

    Args:
        transcripts_file: Input filename

    Returns:
        A pd.DataFrame with the transcript abundances ("id", "count").

    Raises:
        ValueError: When the input file is neither csv or tsv
    """
    cols: list = ["id", "count"]
    if transcripts_file.endswith(".tsv"):
        return pd.read_table(transcripts_file, header=None, names=cols)
    if transcripts_file.endswith(".csv"):
        return pd.read_csv(transcripts_file, header=None, names=cols)
    raise ValueError("File type needs to be either csv or tsv")


def filter_df(gtf_df: pd.DataFrame, transcripts: list) -> pd.DataFrame:
    """Filter dataframe.

    Filter annotations to include only exons
    with the highest transcript support level, i.e. TSL1.

    `feature` column is filtered on value "exon" and
    `free_text` column is filtered to include the string
        denoting the highest transcript support level
    ('transcript_support_level "1"').

    If a list of transcript IDs is given, `free_text` column
        is filtered to include one of the IDs.

    Args:
        df: A pd.DataFrame containing an unparsed gtf-file
        transcript: list of transcript IDs

    Returns:
        A pd.DataFrame containing only rows with exon annotations
            of highest transcript support level and,
        if provided, belonging to one of the given transcripts
    """
    if transcripts is None:
        transcripts = []
    df_filter = gtf_df[
        (gtf_df["feature"] == "exon")
        & (gtf_df["free_text"].str.contains('transcript_support_level "1'))
    ]
    if len(transcripts) > 0:
        df_filter = df_filter[df_filter["free_text"].str.contains(
            "|".join(transcripts), regex=True
        )]

    return df_filter


def str_to_dict(gene_string: str) -> dict:
    """Split between key/value pairs.

    Split string based on delimiter ';' into items, remove empty items and
        split items on delimiter ' ' into
    key/value pairs. Remove quotes from value strings and create a dictionary.

    Args:
        s: A string of the form 'gene_id "GENE1"; transcript_id "TRANSCRIPT1";'

    Returns:
        A dictionary containing e.g. \
            {'gene_id': 'GENE1', 'transcript_id': 'TRANSCRIPT1'}
    """
    # split into items
    # remove empty items
    # split items into key/value pairs
    item_list: list = [x.split() for x in gene_string.split(";") if len(x) > 0]
    # remove quotes for values and return dictionary
    return {item[0]: item[1].strip('"') for item in item_list}


def dict_to_str(gene_dict: dict) -> str:
    """Parse dictionary in gtf free_text column format.

    Takes e.g. dictionary {'gene_id': 'GENE1', 'transcript_id': 'TRANSCRIPT1'}
    and returns string 'gene_id "GENE1"; transcript_id "TRANSCRIPT1";'.
    Key/value pairs are joined by space to form an item and items are
        joinded by ';' to form a string.
    If a value is Not a Number (nan), the key/value pair is omitted
        from the string.

    Args:
        d: A dictionary of the form {'gene_id': 'GENE1', \
            'transcript_id': 'TRANSCRIPT1'}

    Returns:
        A string, e.g. 'gene_id "GENE1"; transcript_id "TRANSCRIPT1";'.
    """
    # join key, value pairs from dictionary with a space in a list,
    # then join items in list by ;
    # end on ;
    # value == value checks that value is not nan
    gene_string: str = "; ".join(
        [f'{key} "{value}"' for key, value in gene_dict.items()]
        ) + ";"
    return gene_string


def reverse_parse_free_text(df_all: pd.DataFrame) -> pd.DataFrame:
    """Reverse parse a gtf based pd.DataFrame.

    The data frame will include only columns that
    are well defnined by gtf-file standards.

    The first 8 defined columns are constant as defined by gtf-file standards.
    Further columns are assumed to be parsed free-text columns
        (see Gtf.parse_free_text()).
    The parsed free-text columns are aggregated as a dictionary and
        the dictionry is parsed as a string in gtf format.

    Args:
        df_all: A pd.DataFrame containing a parsed gtf-file.

    Returns:
        A pd.DataFrame with the columns as defined by gtf-file standards.
    """
    # Define pd.DataFrame containing only parsed free-text columns
    df_free_text = df_all.iloc[:, 8:]
    # Define pd.DataFrame containing only non-parsed columns
    df_non_parsed = df_all.iloc[:, :8]
    # Reverse parsing of free-text columns and add the result as column \
    # `free_text` to output pd.DataFrame
    df_non_parsed["free_text"] = df_free_text.agg(
        pd.Series.to_dict, axis=1
        ).apply(dict_to_str)
    return df_non_parsed


def write_gtf(gtf_df: pd.DataFrame, filename: str) -> None:
    """Save a Gtf object to file in gtf format.

    Makes sure data types are correct and saves object in gtf format.

    Args:
        gtf_df: A pd.DataFrame containing a gtf-file.
        filename: File to save to.
    """
    # Make sure the data types are correct.
    gtf_df = gtf_df.astype(Gtf.dtypes)

    gtf_df.to_csv(
        filename,
        sep="\t",
        header=False,
        index=False,
        quotechar="'",
        mode="a",
    )


class Gtf:
    """Class to read transcripts annotations file into a Gtf object.

    Attributes:
        dtypes: A dictionary containing column names and respective data types.
        parsed: A boolean indicating if the pd.DataFrame is parsed.
        original_columns: A list of columns not touched by parsing.
        free_text_columns: A list of columns created during parsing
            of column `free_text`.
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
        self.data_frame = None
        self.parsed = False
        self.original_columns = list(self.dtypes.keys())
        self.free_text_columns = []

    def read_file(self, annotations_file: str) -> None:
        """Read gtf-file.

        Iterate over chunks of the gtf-file reading 100000 rows at a time.
        Filter chunks for exon annotations of the highest transcript support
        level. Concatenate chunks to get resulting pd.DataFrame.

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
            names=list(self.dtypes.keys()),
            dtype=self.dtypes,
            chunksize=100000,
            iterator=True,
        )
        self.data_frame = pd.concat(
            [filter_df(chunk, transcripts=[]) for chunk in reader]
            )

    def from_dataframe(self, gtf_df: pd.DataFrame) -> None:
        """Initialize Gtf object from pandas Dataframe.

        Part of initialization is:
        Set dataframe attribute
        Check which columns belong to the free-text part of the gtf-file.
        Check if there are no columns called free-text and if so, sets
            the value of parsed attribute to TRUE.

        Args:
            gtf_df: A pd.DataFrame containing a gtf-file.
        """
        self.free_text_columns = [
            col for col in gtf_df.columns if col not in self.original_columns
        ]
        self.data_frame = gtf_df
        if "free_text" not in gtf_df.columns:
            self.parsed = True

    def parse_key_value(self):
        """Parse key/value pairs.

        From `free_text` column into column `key` with row entry `value`.
        Creates a dataframe with columns for keys in the free-text column
        instead of `free_text` column.
        Saves it to Gtf.df attribute.
        """
        assert self.parsed is False
        # create dataframe with columns for values in free_text column
        df_free_text = self.data_frame["free_text"].map(
            str_to_dict
            ).apply(pd.Series)
        # remember which columns come from free_text
        self.free_text_columns = df_free_text.columns
        # join free_text columns to original dataframe and \
        # drop the "free_text" column itself
        self.data_frame = self.data_frame.drop("free_text", axis=1)
        self.original_columns = self.data_frame.columns
        self.data_frame = self.data_frame.join(df_free_text, how="inner")
        # remember that current dataframe is parsed, \
        # i.e. can't be written in gtf format
        self.parsed = True

    def reverse_parse_free_text(self):
        """Reverses parsing of `free_text` column.

        Creates a data frame that can be written in gtf format to file.
        Parsed free-text columns are aggregated
        into `free_text` column according to gtf format specification.
        """
        assert self.parsed is True
        # create dataframe with only free_text columns
        df_free_text = self.data_frame[self.free_text_columns]
        # filter current dataframe to only original columns, \
        # except "free_text" column
        self.data_frame = self.data_frame[self.original_columns]
        # undo parsing and save result in "free_text" column
        self.data_frame["free_text"] = df_free_text.agg(
            pd.Series.to_dict, axis=1
            ).apply(dict_to_str)
        # remember that current dataframe is not parsed
        self.parsed = False

    def pick_transcript(self, transcript_id: str) -> pd.DataFrame:
        """Filter annotations to a given transcript ID."""
        return self.data_frame.query(f"transcript_id == '{transcript_id}'")


class TranscriptGenerator:
    """Class to sample a transcript."""

    # pylint: disable=W0613
    def __new__(
        cls,
        transcript_id: str,
        transcript_count: int,
        transcript_df: pd.DataFrame,
        prob_inclusion: float,
    ):
        """Initialize TranscriptGenerator object."""
        strands = transcript_df["strand"].unique()
        if len(transcript_df) == 0:
            LOG.warning(
                "Transcript \"%s\" can't be sampled: "
                "Annotation is missing or TSL is not 1.", transcript_id
            )
            instance = None
        elif len(strands) > 1:
            LOG.warning(
                "Transcript \"%s\" can't be sampled: Transcript generator is "
                "not implemented for transcripts with exons annotated on "
                "different strands.", transcript_id,
            )
            instance = None
        else:
            instance = super().__new__(cls)

        return instance

    def __init__(
        self,
        transcript_id: str,
        transcript_count: int,
        transcript_df: pd.DataFrame,
        prob_inclusion: float,
    ):
        """Initialize TranscriptGenerator object."""
        self.ts_id = transcript_id
        self.count = transcript_count
        self.data_frame = transcript_df
        self.no_exons = len(transcript_df)
        self.strand = self.data_frame["strand"].unique().item()
        self.prob_inclusion = prob_inclusion

    def get_inclusions(self) -> np.ndarray:
        """Generate inclusions array.

        Each column corresponds to one sample and the number of columns
        corresponds to the number of samples.

        Returns:
            A boolean np.array, where True means intron inclusion.
        """
        inclusion_arr = np.random.rand(
            self.no_exons, self.count
            ) < self.prob_inclusion
        if self.strand == "+":
            inclusion_arr[-1, :] = False
        elif self.strand == "-":
            inclusion_arr[-1, :] = False

        return inclusion_arr

    def get_unique_inclusions(self) -> tuple[list, np.ndarray, np.ndarray]:
        """Get unique inclusions.

        Inclusion of unique intron inclusion via arrays and counts and
        name generation of each unique count.

        Args:

        Returns:
            - List of names for generated exons.
            - A boolean np.array where columns correspond to generated
                transcripts and rows to intron inclusion.
            - A np.array containing sample number per generated inclusions,
                i.e. transcript.
        """
        inclusion_arr = self.get_inclusions()
        # Unique intron inclusion arrays and counts
        inclusion_arr_unique, counts = np.unique(
            inclusion_arr, axis=1, return_counts=True
        )
        # Name for each generated transcript
        names = []
        for i in range(inclusion_arr_unique.shape[1]):
            if np.all(inclusion_arr_unique[:, i] is False, axis=0):
                names.append(self.ts_id)
            else:
                names.append(f"{self.ts_id}_{i}")

        return names, inclusion_arr_unique, counts

    def get_df(
            self, inclusions: np.ndarray, transcript_id: str
            ) -> pd.DataFrame:
        """Get dataframe.

        Take as input a dataframe filtered to one transcript and
        a boolean vector denoting intron inclusions.

        Args:
            inclusions: A boolean vector denoting intron inclusion.
            transcript_id: The transcript id.

        Returns:
            The generated transcript as a pd.DataFrame.
        """
        df_generated = self.data_frame.copy()
        if self.strand == "+":
            original_end = df_generated["end"]
            df_generated["end"] = np.where(
                inclusions,
                df_generated["start"].shift(periods=-1, fill_value=-1) - 1,
                original_end,
            )
        if self.strand == "-":
            original_start = df_generated["start"]
            df_generated["start"] = np.where(
                inclusions,
                df_generated["end"].shift(periods=-1, fill_value=-1) + 1,
                original_start,
            )

        original_id = df_generated["exon_id"]
        df_generated["exon_id"] = np.where(
            inclusions,
            df_generated["exon_id"] + "_" + np.arange(
                len(df_generated)
            ).astype(str),
            original_id,
        )

        df_generated["transcript_id"] = transcript_id
        return df_generated

    def write_sequences(self, filename: str) -> None:
        """Write transcripts to file.

        Args:
            filename: Output csv filename.
        """
        ids, _, counts = self.get_unique_inclusions()
        with open(filename, "a", encoding="utf_8") as file_handle:
            # Add header to output csv for cdna-generator
            if file_handle.tell() == 0:
                file_handle.write(
                    "ID of transcript,ID of parent transcript,"
                    "Transcript copy number\n"
                )

            for transcript_id, transcript_count in zip(ids, counts):
                file_handle.write(
                    f"{transcript_id},{self.ts_id},{transcript_count}\n"
                    )

    def write_annotations(self, filename: str) -> None:
        """Generate a annotations in gtf format for sampled transcript.

        Args:
            Filename: Output gtf-filename.

        Raises:
            ValueError: If given transcript ID could not be sampled.
        """
        ids, inclusions, _ = self.get_unique_inclusions()
        n_unique = len(ids)

        data_frame = pd.concat(
            [self.get_df(inclusions[:, i], ids[i]) for i in range(n_unique)]
        )
        data_frame = reverse_parse_free_text(data_frame)

        write_gtf(data_frame, filename)
        LOG.debug("Transcript \"%s\" sampled.", self.ts_id)


def sample_transcripts(
    input_transcripts_file: str,
    input_annotations_file: str,
    prob_inclusion: float,
    output_transcripts_file: str,
    output_annotations_file: str,
):
    """Sample transcripts.

    Read input files, iterate over transcript IDs,
    sample each transcript and save results.

    Args:
        input_transcripts_file: Filename of transcript abundances,
            needs to be csv or tsv.
        input_annotations_file: Filename of annotations,
            needs to be gtf.
        prob_inclusion: Probability of intron inclusion,
            needs to be float in range [0,1].
        output_transcripts_file: Filename of file to write
            sampled transcripts to.
        output_annotations_file: Filename of file to write
            generated annotations to.
    """
    LOG.info("Probability of intron inclusion: %s", str(prob_inclusion))
    LOG.info("Parsing transcript abundances...")
    transcripts = read_abundances(input_transcripts_file)
    # Make sure that abundance is positive
    transcripts = transcripts[transcripts["count"] > 0]
    LOG.info("Done parsing...")

    LOG.info("Parsing annotations...")
    annotations = Gtf()
    annotations.read_file(input_annotations_file)
    annotations.parse_key_value()
    LOG.info("Done parsing...")

    LOG.info("Start sampling transcripts...")

    for _, row in transcripts.iterrows():
        transcript_id = row["id"]
        transcript_count = row["count"]

        transcript_df = annotations.pick_transcript(transcript_id)
        transcript_generator = TranscriptGenerator(
            transcript_id,
            transcript_count,
            transcript_df,
            prob_inclusion=prob_inclusion,
        )
        try:
            transcript_generator.write_annotations(output_annotations_file)
            transcript_generator.write_sequences(output_transcripts_file)
        except AttributeError:
            pass
    LOG.info("Done.")
