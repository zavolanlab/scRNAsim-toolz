import sys
import warnings
import logging

from cli import parser
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from gtfparse import read_gtf

# ignore warnings from read_gtf
warnings.filterwarnings(action="ignore", category=FutureWarning)


def compliment(res: str) -> str:
    translate_dict = {"A": "T", "T": "A", "U": "A", "G": "C", "C": "G"}
    if res not in translate_dict.keys():
        print(f"Unknown character, {res}")
        sys.exit(1)
    return translate_dict[res]


def seq_compliment(sequence: str) -> str:
    if sequence is None:
        return "None"
    _ = "".join([compliment(char) for char in str(sequence)])[::-1]  # reverse string
    return _


class CDNAGen:
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

        self.run()

    def run(self):
        self.read_csv()
        self.read_fasta()
        self.read_gtf()
        self.add_sequences()
        self.add_compliment()
        self.add_records()
        self.write_fasta()
        self.write_csv()

    def add_records(self):
        self.fasta_records = []
        for index, row in self.df_input_GTF.iterrows():
            if row["compliment"] is not None:
                copy_number = row["Transcript_Copy_Number"]
                record = SeqRecord(
                    Seq(row["compliment"]),
                    row["cdna_ID"],
                    f"Transcript copy number: {copy_number}",
                    "",
                )
                self.fasta_records.append(record)

    def add_sequences(self):
        self.df_input_GTF["priming_site"] = self.df_input_GTF.apply(
            lambda row: self.read_primingsite(row["seqname"], row["start"]),
            axis=1,
        )

    def add_compliment(self):
        self.df_input_GTF["compliment"] = self.df_input_GTF["priming_site"].apply(
            lambda x: seq_compliment(x)
        )

    def read_primingsite(self, sequence, start):
        if sequence not in self.fasta_dict.keys():
            return None
        _ = self.fasta_dict[sequence].seq[start:]
        return _

    def read_fasta(self):
        record = SeqIO.parse(self.fasta, "fasta")
        records = list(record)
        self.fasta_dict = {x.name: x for x in records}

    def read_csv(self):
        df_input_CSV = pd.read_csv(self.cpn, index_col=False)
        df_input_CSV = (
            df_input_CSV.reset_index()
        )  # make sure indexes pair with number of rows
        self.df_input_CSV = df_input_CSV

    def read_gtf(self):
        # returns GTF with essential columns such as "feature", "seqname", "start", "end"
        # alongside the names of any optional keys which appeared in the attribute column
        df_input_GTF = read_gtf(self.gtf)
        df_input_GTF["Binding_Probability"] = pd.to_numeric(
            df_input_GTF["Binding_Probability"]
        )  # convert to numeric
        df_normalization_bind_probablility = df_input_GTF.groupby("seqname")[
            "Binding_Probability"
        ].sum()  # extract binding probablility
        count = 0
        prev_id = None
        # Adds Normalized_Binding_Probability and Transcript_Copy_Number to each transcript in the dataframe
        for index, row in df_input_GTF.iterrows():
            # GTF transcript ID
            id_GTF = str(row["seqname"])
            if id_GTF == prev_id:
                count += 1
            else:
                prev_id = None
                count = 0
                # CVS transcript ID
            id_CSV = str(row["seqname"]).split("_")[1]
            # Calculate Normalized_Binding_Probability and add to GTF dataframe
            df_input_GTF.loc[index, "Normalized_Binding_Probability"] = (
                row["Binding_Probability"] / df_normalization_bind_probablility[id_GTF]
            )
            # Calculate Normalized_Binding_Probability and add to GTF dataframe
            csv_transcript_copy_number = self.df_input_CSV.loc[
                self.df_input_CSV["ID of transcript"] == int(id_CSV),
                "Transcript copy number",
            ].iloc[0]
            df_input_GTF.loc[index, "Transcript_Copy_Number"] = round(
                csv_transcript_copy_number
                * df_input_GTF.loc[index, "Normalized_Binding_Probability"]
            )
            df_input_GTF.loc[index, "cdna_ID"] = f"{id_GTF}_{count}"
            prev_id = id_GTF

        self.df_input_GTF = df_input_GTF

    def write_fasta(self):
        print(self.fasta_records)
        SeqIO.write(self.fasta_records, self.output_fasta, "fasta")

    def write_csv(self):
        self.df_input_GTF[["cdna_ID", "Transcript_Copy_Number"]].to_csv(
            self.output_csv, index=False
        )

    def return_output(self):
        return self.output_fasta, self.output_csv


if __name__ == "main":
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO,
    )
    LOG = logging.getLogger(__name__)
    cnda_object = parser()
