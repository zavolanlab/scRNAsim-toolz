import sys
import warnings

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from gtfparse import read_gtf

# ignore warnings from read_gtf
warnings.filterwarnings(action="ignore", category=FutureWarning)


def compliment(res: str) -> str:
    """
    Returns the cDNA compliment of a given base pair
    Args:
        res: residue code.

    Returns: corresponding cDNA residue.

    """
    translate_dict = {"A": "T", "T": "A", "U": "A", "G": "C", "C": "G"}
    if res not in translate_dict.keys():
        print(f"Unknown character, {res}")
        sys.exit(1)
    return translate_dict[res]


def seq_compliment(sequence: str) -> str:
    """
    Returns the corresponding cDNA sequence by finding the complimentary
    base pairs and returning the reversed sequence.

    Args:
        sequence: sequence to be converted into cDNA.

    Returns: corresponding cDNA sequence.

    """
    if sequence is None:
        return "None"
    _ = "".join([compliment(char) for char in str(sequence)])[::-1]  # reverse string
    return _


class CDNAGen:
    """
    Module that performs the cDNA synthesis.
    """
    def __init__(self, ifasta: str, igtf: str, icpn: str, ofasta: str, ocsv: str):
        # inputs
        self.fasta = ifasta
        self.gtf = igtf
        self.cpn = icpn
        self.output_fasta = ofasta
        self.output_csv = ocsv

        # variables
        self.fasta_dict = None
        self.fasta_records = None
        self.gtf_df = None
        self.run()

    def run(self) -> None:
        """
        Executes the cDNA workflow.
        Returns: None

        """
        self.read_csv()
        self.read_fasta()
        self.read_gtf()
        self.add_sequences()
        self.add_compliment()
        self.add_records()
        print()  # blank line for pretty printing
        self.write_fasta()
        self.write_csv()

    def add_records(self) -> None:
        self.fasta_records = []
        for index, row in self.gtf_df.iterrows():
            if row["compliment"] is not None:
                copy_number = row["Transcript_Copy_Number"]
                record = SeqRecord(
                    Seq(row["compliment"]),
                    row["cdna_ID"],
                    f"Transcript copy number: {copy_number}",
                    "",
                )
                self.fasta_records.append(record)

    def add_sequences(self) -> None:
        """
        Adds the sequence for a given priming site.
        Returns: None

        """
        self.gtf_df["priming_site"] = self.gtf_df.apply(
            lambda row: self.read_primingsite(row["seqname"], row["start"]),
            axis=1,
        )

    def add_compliment(self) -> None:
        """
        Adds the complimentary cDNA sequence.
        Returns: None

        """
        self.gtf_df["compliment"] = self.gtf_df["priming_site"].apply(
            lambda x: seq_compliment(x)
        )

    def read_primingsite(self, sequence: str, start: int) -> None:
        """Read a fasta file from a given start character

        Reads a fasta sequence with ID (sequence) and returns the
        sequence starting from the index start.

        Args:
            sequence: sequence ID to be read.
            start: start of the sequence.

        Returns: None

        """
        if sequence not in self.fasta_dict.keys():
            return None
        _ = self.fasta_dict[sequence].seq[start:]
        return _

    def read_fasta(self) -> None:
        """Read a given fasta file.

        Wrapper for SeqIO.parse.

        Returns: None

        """
        record = SeqIO.parse(self.fasta, "fasta")
        records = list(record)
        self.fasta_dict = {x.name: x for x in records}

    def read_csv(self) -> None:
        """Reads a given copy number csv file

        Wrapper for Pandas read_csv.

        Returns: None

        """
        df_csv = pd.read_csv(self.cpn, index_col=False)
        df_csv = (
            df_csv.reset_index()
        )  # make sure indexes pair with number of rows
        self.csv_df = df_csv

    def read_gtf(self) -> None:
        """Read and process the GTF file.

        Reads a GTF file and determines copy numbers from normalized probabilities.

        Returns: None

        """
        # returns GTF with essential columns such as "feature", "seqname", "start", "end"
        # alongside the names of any optional keys which appeared in the attribute column
        gtf_df = read_gtf(self.gtf)
        gtf_df["Binding_Probability"] = pd.to_numeric(
            gtf_df["Binding_Probability"]
        )  # convert to numeric
        df_normalization_bind_probablility = gtf_df.groupby("seqname")[
            "Binding_Probability"
        ].sum()  # extract binding probability
        count = 0
        prev_id = None
        # Adds Normalized_Binding_Probability and Transcript_Copy_Number to each transcript in the dataframe
        for index, row in gtf_df.iterrows():
            # GTF transcript ID
            id_ = str(row["seqname"])
            if id_ == prev_id:
                count += 1
            else:
                prev_id = None
                count = 0
                # CVS transcript ID
            id_csv = str(row["seqname"]).split("_")[1]
            # Calculate Normalized_Binding_Probability and add to GTF dataframe
            gtf_df.loc[index, "Normalized_Binding_Probability"] = (
                row["Binding_Probability"] / df_normalization_bind_probablility[id_]
            )
            # Calculate Normalized_Binding_Probability and add to GTF dataframe
            csv_transcript_copy_number = self.csv_df.loc[
                self.csv_df["ID of transcript"] == int(id_csv),
                "Transcript copy number",
            ].iloc[0]
            gtf_df.loc[index, "Transcript_Copy_Number"] = round(
                csv_transcript_copy_number
                * gtf_df.loc[index, "Normalized_Binding_Probability"]
            )
            gtf_df.loc[index, "cdna_ID"] = f"{id_}_{count}"
            prev_id = id_

        self.gtf_df = gtf_df

    def write_fasta(self) -> None:
        """Writes cDNA fasta records to file.

        Wrapper for SeqIO.write.

        Returns: None

        """
        SeqIO.write(self.fasta_records, self.output_fasta, "fasta")
        print(f"Fasta file successfully written to: {self.output_fasta}")

    def write_csv(self) -> None:
        """Writes the copy number information to a csv file.

        Wrapper for Pandas to_csv.

        Returns: None

        """
        self.gtf_df[["cdna_ID", "Transcript_Copy_Number"]].to_csv(
            self.output_csv, index=False
        )
        print(f"Copy number csv file successfully written to: {self.output_csv}")
