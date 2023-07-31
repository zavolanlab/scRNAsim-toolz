"""Find representative transcripts."""

import logging

LOG = logging.getLogger(__name__)


class FindRepTrans:
    """Find representative transcripts."""

    def __init__(self):
        """Initiate."""

    @staticmethod
    def attributes_converter(attributes: str) -> list:
        """Attributes converter function.

        This funtion converts the "unstructured" ;-seperated part of
        the line into a list of identifiers and corresponding data,
        the structure of which can be used ot find the data easily e.g
        the index of the identifier transcript_id + 1 will give the
        transcript id of the current gene.
        Input:
            attributes = str() # the unstructured part of the entry
        Output:
            attributes = list() # cleaned list with the \
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
        else:
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
        # setting default variables
        rep_transcripts = dict()
        cur_g_id = ""
        # [transcript_id, transcript_support_level, transcript_length]
        cur_best_trans = ["", 100, 0]

        with open(file_name, "r", encoding="utf-8") as file:
            for line in file:
                entry = line.split("\t")

                # removes expected but unneeded entries
                if len(entry) == 1 or entry[2] in [
                    "CDS", "stop_codon",
                    "five_prime_utr", "three_prime_utr",
                    "start_codon", "Selenocysteine"
                        ]:
                    continue

                # this function turns the less organized part of the entry
                # into a readable list
                attributes = self.attributes_converter(entry[8])

                # looking for and processing exons entries
                if entry[2] == "exon":
                    if cur_g_id != attributes[1]:
                        LOG.error("Exon from an unexpected gene")
                        raise ValueError("Exon from an unexpected gene")
                    elif (
                        self.find_in_attributes(
                            attributes, "transcript_id"
                        ) != cur_tID
                        ):
                        LOG.error("Exon from an unexpected transcript")
                        raise ValueError("Exon from an unexpected transcript")

                    # adding the length of the exon to the appropriate list and
                    # checking for changes in best transcript
                    if pot_best_trans:
                        pot_best_trans[2] += int(entry[4]) - int(entry[3])
                        if pot_best_trans[2] > cur_best_trans[2]:
                            cur_best_trans = pot_best_trans
                            pot_best_trans = False
                    else:
                        cur_best_trans[2] += int(entry[4]) - int(entry[3])

                # looking for and processing transcript entries
                elif entry[2] == "transcript":
                    # verify that the gen is correct
                    if cur_g_id != attributes[1]:
                        LOG.error("Transcript from an unexpected gene")
                        raise ValueError("Transcript from an unexpected gene")

                    # finding the transcript id and the support level
                    cur_tID = self.find_in_attributes(
                        attributes, "transcript_id"
                        )
                    t_supp_lvl = self.find_in_attributes(
                        attributes, "transcript_support_level"
                        )

                    # If there is no transcript support level or the level is
                    # given as NA it is nomed as 100. else the transcript
                    # support level is turned into int
                    if t_supp_lvl == "NA":
                        t_supp_lvl = 100
                    else:
                        if t_supp_lvl.isdigit():
                            t_supp_lvl = int(t_supp_lvl)
                        else:
                            t_supp_lvl = 100

                    # decides if the transcript has potential to become the
                    # representative transcript
                    if t_supp_lvl < cur_best_trans[1] or cur_best_trans[0] == "":
                        cur_best_trans = [cur_tID, t_supp_lvl, 0]
                        pot_best_trans = False
                        ignor_trans = False
                    elif t_supp_lvl == cur_best_trans[1]:
                        pot_best_trans = [cur_tID, t_supp_lvl, 0]
                    else:
                        ignor_trans = True

                # looking for and processing gene entries
                elif entry[2] == "gene":
                    # updating rep_transcripts dict
                    if cur_g_id in rep_transcripts:
                        if (rep_transcripts[cur_g_id][1] > cur_best_trans[1]
                            or (rep_transcripts[cur_g_id][1] ==
                                cur_best_trans[1]
                                and rep_transcripts[cur_g_id][2] <
                                cur_best_trans[2])):
                            rep_transcripts[cur_g_id] = cur_best_trans
                    else:
                        rep_transcripts[cur_g_id] = cur_best_trans

                    # updating cur_g_id and resetting cur_best_trans
                    cur_g_id = attributes[1]
                    cur_best_trans = ["", 100, 0]

                # raises an error for unidentifiable entries
                else:
                    LOG.error("This entry could not be identified")
                    raise ValueError("This entry could not be identified")

            # adding the final gene to the dictionary
            if cur_g_id in rep_transcripts:
                if (rep_transcripts[cur_g_id][1] > cur_best_trans[1] or
                        (rep_transcripts[cur_g_id][1] == cur_best_trans[1] and
                        rep_transcripts[cur_g_id][2] < cur_best_trans[2])):
                    rep_transcripts[cur_g_id] = cur_best_trans
            else:
                rep_transcripts[cur_g_id] = cur_best_trans

            del rep_transcripts[""]
            rep_transcripts = self.reformat_reptrans(rep_transcripts)
            return rep_transcripts

    def gtf_file_writer(self, original_file: str,
                        rep_transcript_dict: dict, output_file: str):
        """Gtf file writer.

        This function writes the output GTF file.
        """
        output = []

        with open(original_file, "r", encoding="utf-8") as f:
            for line in f:
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


# def _test():
#     """
#     This funtion is meant to be run for test
#     Output:
#         file with the dictionary generated based on the test file
#     """
#     file_name = "test.gtf"
#     rt = get_rep_trans(file_name)
#     expected_result = {"ENSG00000160072": "ENST00000472194",
#                        "ENSG00000234396": "ENST00000442483",
#                        "ENSG00000225972": "ENST00000416931",
#                        "ENSG00000224315": "ENST00000428803",
#                        "ENSG00000198744": "ENST00000416718",
#                        "ENSG00000279928": "ENST00000624431",
#                        "ENSG00000228037": "ENST00000424215",
#                        'ENSG00000142611': 'ENST00000378391'}
#     if rt != expected_result:
#         print("The test failed due to not yielding the same results")
#         print("The results the program got\n", rt)
#         print("The expected results\n", expected_result)
#     else:
#         print("The test was successful")


# # Execution part #
# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(
#         description="find_representativ_transcripts",
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter
#         )
#     parser.add_argument("-file_name", required=True,
#                         help="gtf file with genome annotation")
#     parser.add_argument("-t", required=False, default=False,
#                         help="to run the test input -t True")
#     args = parser.parse_args()

#     # standadize the file_name inlude .gtf#
#     file_name = args.file_name
#     i_gtf = file_name.find(".gtf")
#     if i_gtf == -1:
#         file_name += ".gtf"

#     if args.t:
#         _test()
#     else:
#         get_rep_trans(file_name)
