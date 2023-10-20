"""Find representative transcripts."""
import logging
from typing import Union
import pandas as pd  # type: ignore

LOG = logging.getLogger(__name__)


# pylint: disable=R0912,R0915
class FindRepTrans:
    """Find representative transcripts."""

    def __init__(self):
        """Initiate."""

    @staticmethod
    def attributes_converter(attributes):
        """Attributes converter function.

        This funtion converts the "unstructured" ;-seperated part of
        the line into a list of identifiers and corresponding data,
        the structure of which can be used ot find the data easily e.g
        the index of the identifier transcript_id + 1 will give the
        transcript id of the current gene.
        Input:
            attributes = str() # the unstructured part of the entry
        Output:
            attributes = list() # cleaned list with the
                                  characteristics described above
        """
        attributes = (
            attributes.replace('"', "")
            .replace(";", "")
            .replace("\\n", "")
            .split(" ")
        )
        return attributes

    @staticmethod
    def find_in_attributes(attributes: list, look_for: str) -> str:
        """Find in attributes function.

        This function finds a keyword and used that to locate the value of that
        keyword e.g key = gene_id, value = 'ENSMUSG00002074970',
        this works as they are next to each other in the attributes list.
        Inputs:
            attributes = list()
            look_for = str() # string of the name of the key to look for
        Output:
            attributes[index] or NA = str() # NA is returned if the key
                                            was not found in the attributes
        """
        if look_for in attributes:
            index = attributes.index(look_for) + 1
            return attributes[index]
        LOG.warning('No %s in the entry, the return was set to NA',
                    look_for)
        return "NA"

    @staticmethod
    def reformat_reptrans(rep_trans_dict: dict) -> dict:
        """Reformat dictionary.

        This function is meant to reformat dictionary of the representative
        transcripts into an dictionary with only one entry per key
        Input:
            rep_trans_dict = {gene_id : [
                transcript_id, transcript_support_level, transcript_length]}
        Output:
            rep_transcripts = {gene_id : transcript_id}
        """
        rep_transcripts = {}
        for gene_id in rep_trans_dict:
            rep_transcripts[gene_id] = rep_trans_dict[gene_id][0]

        return rep_transcripts

    def get_rep_trans(self, file_name: str) -> dict:
        """Get representative transcripts.

        This is the main function of this script. It selects one
        representative transcript per gene based on a GTF annotation file.
        It does so by two criteria: the transcript support level and if
        there are several transcripts of one gene that have the same
        transcript_support_level, it chooses the one that corresponds
        to the longest mRNA.

        Args:
            file_name (str): Name of the annotation file with or without
            the .gtf extension.

        Returns:
            rep_transcripts (dict): Dictionary of gene_id to transcript_id
            representing the selected representative transcripts.

        Raises:
            ValueError: If an unexpected entry is encountered in the GTF file.
        """
        result_list = []

        with open(file_name, "r", encoding="utf-8") as file:
            for line in file:
                entry = line.strip().split("\t")

                if len(entry) == 1 or entry[2] in [
                    "CDS", "stop_codon",
                    "five_prime_utr", "three_prime_utr",
                    "start_codon", "Selenocysteine"
                ]:
                    continue

                attributes = self.attributes_converter(entry[8])

                # 1. GENE entries
                if entry[2] == "gene":
                    cur_g_id = self.find_in_attributes(attributes, "gene_id")

                # 2. TRANSCRIPT entries
                elif entry[2] == "transcript":
                    cur_t_id = self.find_in_attributes(
                        attributes, "transcript_id"
                        )
                    tsl: Union[int, str] = self.find_in_attributes(
                        attributes, "transcript_support_level"
                        )
                    if tsl == "NA":
                        tsl = 100
                    else:
                        tsl = int(tsl)

                    result_list.append(
                        {"gene_id": cur_g_id, "transcript_id": cur_t_id,
                         "tsl": tsl, "exon_length": 0}
                    )

                # EXON entries
                elif entry[2] == "exon":
                    exon_length = int(entry[4]) - int(entry[3])

                    for item in result_list:
                        if (
                            item["gene_id"] == cur_g_id
                            and item["transcript_id"] == cur_t_id
                        ):
                            item["exon_length"] += exon_length  # type: ignore

        # Convert the list of dictionaries into a DataFrame
        df = pd.DataFrame(result_list)

        # Create dictionary of representative transcripts
        sorted_df = df.sort_values(
            by=['gene_id', 'tsl', 'exon_length'], ascending=[True, True, False]
            )

        # Keep the first row for each gene_id (best representative transcript)
        best_rep_df = sorted_df.drop_duplicates(subset='gene_id', keep='first')

        # Create a dictionary with gene_id as keys and transcript_id as values
        best_rep_transcripts = best_rep_df.set_index(
            'gene_id'
        )['transcript_id'].to_dict()

        return best_rep_transcripts

    def gtf_file_writer(self, original_file: str,
                        rep_transcript_dict: dict, output_file: str):
        """Gtf file writer.

        This function writes the output GTF file.
        """
        output = []

        with open(original_file, "r", encoding="utf-8") as file:
            for line in file:
                if line.startswith("#"):
                    continue

                entry = line.split("\t")
                attributes = self.attributes_converter(entry[8])
                feature_type = entry[2]

                if feature_type == "gene":
                    gene_id = self.find_in_attributes(attributes, "gene_id")
                    output.append(line)
                else:
                    transcript_id = self.find_in_attributes(
                        attributes, "transcript_id"
                        )
                    if gene_id in rep_transcript_dict and \
                            rep_transcript_dict[gene_id] == transcript_id:
                        output.append(line)

        with open(output_file, "w", encoding="utf-8") as last_file:
            last_file.writelines(output)
