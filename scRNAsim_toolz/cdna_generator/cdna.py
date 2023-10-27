"""cDNA generator."""
import warnings
import logging
from typing import Optional, List, Dict, Any
import pandas as pd  # type: ignore
from Bio import SeqIO  # type: ignore
from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from gtfparse import read_gtf  # type: ignore

LOG = logging.getLogger(__name__)

# ignore warnings from read_gtf
warnings.filterwarnings(action="ignore", category=FutureWarning)


def complement(res: str) -> str:
    """Return the cDNA complement of a given base pair.

    Args:
        res: residue code.

    Returns: corresponding cDNA residue.
    Raises: Value error

    """
    translate_dict = {"A": "T", "T": "A", "U": "A", "G": "C", "C": "G"}
    if res not in translate_dict:
        LOG.warning("Unknown character, %s", res)
        raise ValueError
    return translate_dict[res]


def seq_complement(sequence: str) -> Optional[str]:
    """Return the corresponding cDNA sequence.

    Find the complementary base pairs and
    returning the reversed sequence.

    Args:
        sequence: sequence to be converted into cDNA.

    Returns: corresponding cDNA sequence.

    """
    if sequence is None:
        return None
    _ = "".join([
        complement(char) for char in str(sequence)
        ])[::-1]  # reverse string
    return _


# pylint: disable=R0902
class CDNAGen:
    """Perform the cDNA synthesis."""

    # pylint: disable=R0913
    def __init__(
        self, ifasta: str, igtf: str, icpn: str, ofasta: str, ocsv: str
    ):
        """Initialise function."""
        # inputs
        self.fasta = ifasta
        self.gtf = igtf
        self.cpn = icpn
        self.output_fasta = ofasta
        self.output_csv = ocsv

        # variables
        self.csv_df = pd.DataFrame()
        self.fasta_dict: Dict[str, Any] = {}
        self.fasta_records: List[SeqRecord] = []
        self.gtf_df = pd.DataFrame()
        self.run()

    def run(self) -> None:
        """Execute the cDNA workflow.

        Returns: None

        """
        self.process_csv()
        self.process_fasta()
        self.process_gtf()
        self.add_sequences()
        self.add_complement()
        self.add_records()
        self.write_fasta()
        self.write_csv()

    def process_csv(self) -> None:
        """Read a given copy number csv file.

        Wrapper for Pandas read_csv.

        Returns: None

        """
        df_csv = pd.read_csv(self.cpn, index_col=False)
        df_csv = df_csv.reset_index()  # make sure indexes pair with number of rows # noqa: E501
        self.csv_df = df_csv

    def process_fasta(self) -> None:
        """Read a given fasta file.

        Wrapper for SeqIO.parse.

        Returns: None

        """
        records = list(SeqIO.parse(self.fasta, "fasta"))
        self.fasta_dict = {x.name: x for x in records}

    def process_gtf(self) -> None:
        """Read and process the GTF file.

        Reads a GTF file and determines copy numbers from
        normalized probabilities.

        Returns: None

        """
        # returns GTF with essential columns such as
        # "feature", "seqname", "start", "end"
        # alongside the names of any optional keys
        # which appeared in the attribute column
        gtf_df = read_gtf(self.gtf, result_type="pandas")  # from gtfparse

        gtf_df["Binding_Probability"] = pd.to_numeric(
            gtf_df["Binding_Probability"]
        )  # convert to numeric
        df_norm_bind_prob = gtf_df.groupby("seqname")[
            "Binding_Probability"
        ].sum()  # extract binding probability
        count = 0
        prev_id = None
        # Adds Normalized_Binding_Probability and Transcript_Copy_Number
        # to each transcript in the dataframe
        for index, row in gtf_df.iterrows():
            # GTF transcript ID
            id_ = str(row["seqname"])
            if id_ == prev_id:
                count += 1
            else:
                count = 0  # reset count
            # CSV transcript ID
            id_csv = f"{id_}_{count}"
            # Calculate Normalized_Binding_Probability and add to GTF dataframe
            gtf_df.loc[index, "Normalized_Binding_Probability"] = (
                row["Binding_Probability"] / df_norm_bind_prob[id_]
            )
            # Calculate Normalized_Binding_Probability and add to GTF dataframe
            csv_transcript_copy_number = self.csv_df.loc[
                self.csv_df["ID of transcript"] == id_csv,
                "Transcript copy number",
            ].iloc[0]  # pop the first value in the frame
            gtf_df.loc[index, "Transcript_Copy_Number"] = round(
                csv_transcript_copy_number
                * gtf_df.loc[index, "Normalized_Binding_Probability"]
            )
            gtf_df.loc[index, "cdna_ID"] = f"{id_}_{count}"
            prev_id = id_

        gtf_df['Transcript_Copy_Number'] = gtf_df[
            'Transcript_Copy_Number'
        ].astype(int)
        self.gtf_df = gtf_df

    def add_sequences(self) -> None:
        """Add the sequence for a given priming site.

        Returns: None

        """
        self.gtf_df["priming_site"] = self.gtf_df.apply(
            lambda row: self.read_primingsite(row["seqname"], row["start"]),
            axis=1,
        )

    def read_primingsite(self, sequence: str, end: int) -> None:
        """Read a fasta file from a given start character.

        Reads a fasta sequence with ID (sequence) and returns the
        sequence starting from the index start.

        Args:
            sequence: sequence ID to be read.
            end: end index of the priming site.

        Returns: None

        """
        if sequence not in self.fasta_dict.keys():
            return None
        return self.fasta_dict[sequence].seq[:end]

    def add_complement(self) -> None:
        """Add the complementary cDNA sequence.

        Returns: None

        """
        self.gtf_df["complement"] = self.gtf_df["priming_site"].apply(
            seq_complement
            )

    def add_records(self) -> None:
        """Add data records to fasta file.

        Adds the copy number information to the fasta records.

        Returns: None

        """
        self.fasta_records = []
        for _, row in self.gtf_df.iterrows():
            if row["complement"] is not None:
                copy_number = row["Transcript_Copy_Number"]
                for _ in range(int(copy_number)):
                    record = SeqRecord(
                        Seq(row["complement"]),
                        row["cdna_ID"],
                        f"Transcript copy number: {copy_number}",
                        "",
                    )
                    self.fasta_records.append(record)

    def write_fasta(self) -> None:
        """Write cDNA fasta records to file.

        Wrapper for SeqIO.write.

        Returns: None

        """
        SeqIO.write(self.fasta_records, self.output_fasta, "fasta")
        LOG.info("Fasta file successfully written to: %s", self.output_fasta)

    def write_csv(self) -> None:
        """Write the copy number information to a csv file.

        Wrapper for Pandas to_csv.

        Returns: None

        """
        df_to_save = self.gtf_df[["cdna_ID", "Transcript_Copy_Number"]]
        # Stop outputting header
        df_to_save.to_csv(self.output_csv, index=False, header=False)
        LOG.info("Copy number csv file successfully written to: %s",
                 self.output_csv)
