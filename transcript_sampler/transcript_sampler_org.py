import pandas as pd
import numpy as np
import logging
from gtfparse import read_gtf

LOG = logging.getLogger(__name__)

def attributes_converter(attributes: str) -> list:
    """
    This funtion converts the "unstructured" ;-seperated part of he line into
    a list of identifiers and corresponding data, the structure of
    which can be used ot find the data easily e.g the index of the identifier
    transcript_id + 1 will give the transcript id of the current gene
    Input:
        attributes = str() # the unstructured part of the entry
    Output:
        attributes = list() # cleaned list with the characteristics described
    """
    attributes = (
        attributes.replace('"', "")
        .replace(";", "")
        .replace("\\n", "")
        .split(" ")
    )
    return attributes


def find_in_attributes(attributes: list, look_for: str) -> str:
    """
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
        LOG.warning(f'No {look_for} in the entry, the return was set to NA')
        return "NA"


def _re_format(rep_trans_dict: dict) -> dict:
    """
    This function is meant to reformat dictionary of the representative
    transcripts into an dictionary with only one entry per key
    Input:
        rep_trans_dict = {gene_id : [
            transcript_id, transcript_support_level, transcript_length]}
    Output:
        rep_transcripts = {gene_id : transcript_id}
    """
    rep_transcripts = dict()
    for gene_id in rep_trans_dict:
        rep_transcripts[gene_id] = rep_trans_dict[gene_id][0]

    return rep_transcripts


def get_rep_trans(file_name: str = "test.gtf") -> dict:
    """
    This is the main function of this script. It selects one representative transcript per gene based on a GTF annotation file.
    It does so by two criteria: the transcript support level and if there are several transcripts of one gene that have the same transcript_support_level, it chooses the one that corresponds to the longest mRNA.

    Args:
        file_name (str): Name of the annotation file with or without the .gtf extension.

    Returns:
        rep_transcripts (dict): Dictionary of gene_id to transcript_id representing the selected representative transcripts.

    Raises:
        ValueError: If an unexpected entry is encountered in the GTF file.
    """

    # setting default variables
    rep_transcripts = {}
    cur_gID = ""
    cur_best_trans = ["", 100, 0]  # [transcript_id, transcript_support_level, transcript_length]

    with open(file_name, "r") as f:
        for line in f:
            entry = line.split("\t")

            # removes expected but unneeded entries
            if len(entry) == 1 or entry[2] in [
                "CDS",
                "stop_codon",
                "five_prime_utr",
                "three_prime_utr",
                "start_codon",
                "Selenocysteine"
                ]:
                continue

            # this function turns the less organized part of the entry
            # into a readable list
            attributes = attributes_converter(entry[8])

            # looking for and processing exons entries
            if entry[2] == "exon":
                if ignor_trans:
                    continue
                elif cur_gID != attributes[1]:
                    LOG.error()
                    raise ValueError("Exon from an unexpected gene")
                elif find_in_attributes(attributes, "transcript_id") != cur_tID:
                    LOG.error()
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
                if cur_gID != attributes[1]:
                    LOG.error()
                    raise ValueError("Transcript from an unexpected gene")

                # finding the transcript id and the support level
                cur_tID = find_in_attributes(attributes, "transcript_id")
                t_supp_lvl = find_in_attributes(attributes, "transcript_support_level")

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
                if cur_gID in rep_transcripts:
                    if rep_transcripts[cur_gID][1] > cur_best_trans[1] or (rep_transcripts[cur_gID][1] == cur_best_trans[1] and rep_transcripts[cur_gID][2] < cur_best_trans[2]):
                        rep_transcripts[cur_gID] = cur_best_trans
                else:
                    rep_transcripts[cur_gID] = cur_best_trans

                # updating cur_gID and resetting cur_best_trans
                cur_gID = attributes[1]
                cur_best_trans = ["", 100, 0]

            # raises an error for unidentifiable entries
            else:
                LOG.error()
                raise ValueError("This entry could not be identified")

         # adding the final gene to the dictionary
        if cur_gID in rep_transcripts:
            if rep_transcripts[cur_gID][1] > cur_best_trans[1] or (rep_transcripts[cur_gID][1] == cur_best_trans[1] and rep_transcripts[cur_gID][2] < cur_best_trans[2]):
                rep_transcripts[cur_gID] = cur_best_trans
        else:
            rep_transcripts[cur_gID] = cur_best_trans

        del rep_transcripts[""]
        rep_transcripts = _re_format(rep_transcripts)
        return rep_transcripts


def _test():
    """
    This funtion is meant to be run for test
    Output:
        file with the dictionary generated based on the test file
    """
    file_name = "test.gtf"
    rt = get_rep_trans(file_name)
    expected_result = {
        "ENSG00000160072": "ENST00000472194",
        "ENSG00000234396": "ENST00000442483",
        "ENSG00000225972": "ENST00000416931",
        "ENSG00000224315": "ENST00000428803",
        "ENSG00000198744": "ENST00000416718",
        "ENSG00000279928": "ENST00000624431",
        "ENSG00000228037": "ENST00000424215",
        "ENSG00000142611": "ENST00000378391",
    }
    if rt != expected_result:
        print("The test fail due to not yieding the same results")
        print("The results the program got\n", rt)
        print("The expected results\n", expected_result)
    else:
        print("The test was succsesfull")


def gtf_file_writer(original_file: str, rep_transcript_dict: dict, output_file: str):
    """
    This function writes the output GTF file.
    """
    output = []

    with open(original_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            entry = line.split("\t")
            attributes = attributes_converter(entry[8])
            feature_type = entry[2]

            if feature_type == "gene":
                gene_id = find_in_attributes(attributes, "gene_id")
                output.append(line)
            else:
                transcript_id = find_in_attributes(attributes, "transcript_id")
                if gene_id in rep_transcript_dict and rep_transcript_dict[gene_id] == transcript_id:
                    output.append(line)

    with open(output_file, "w") as last_file:
        last_file.writelines(output)


def gtf_to_df(gtf_file: str) -> pd.DataFrame:
    """
    This function takes a .gtf file and converts it into a pandas DataFrame
    containing gene_id and their transcript_id.
    
    Args:
        gtf_file (str): Path to the .gtf file.
    
    Returns:
        df_gtf (pd.DataFrame): Pandas DataFrame containing columns 'Gene' and 'Transcript'.
    
    Raises:
        None
    """
    df_gtf = read_gtf(gtf_file,).to_pandas()
    df_gtf = df_gtf[df_gtf["feature"] == "transcript"]
    df_gtf = df_gtf[["gene_id", "transcript_id"]]
    df_gtf = df_gtf.rename(columns={"gene_id": "Gene", "transcript_id": "Transcript"})
    return df_gtf


def dict_reprTrans_to_df(dict_reprTrans: "dict[str, str]") -> pd.DataFrame:
    """
    Convert a dictionary of genes and their representative transcript into a DataFrame.
    
    Args:
        dict_reprTrans (dict): {'Gene': ['transcriptA', 'transcriptB'], ...}
    
    Returns:
        Pandas DataFrame with 'Gene' and 'Transcript' as columns.
    
    Raises:
        TypeError: Only dictionaries are allowed.
        TypeError: Keys should be strings.
        TypeError: Values should be strings.
    """
    if not isinstance(dict_reprTrans, dict):
        LOG.error()
        raise TypeError("Only dictionaries are allowed")
    if not all(isinstance(key, str) for key in dict_reprTrans.keys()):
        LOG.error()
        raise TypeError("Keys should be strings")
    if not all(isinstance(value, str) for value in dict_reprTrans.values()):
        LOG.error()
        raise TypeError("Values should be strings")

    df_reprTrans = pd.DataFrame.from_dict(dict_reprTrans, orient="index", columns=["reprTranscript"])
    df_reprTrans = df_reprTrans.reset_index()
    df_reprTrans.columns = ["Gene", "reprTrans"]
    df_reprTrans["reprTrans"] = df_reprTrans["reprTrans"].str.replace(r"\.[1-9]", "", regex=True)

    return df_reprTrans


def tsv_or_csv_to_df(input_txt: str) -> pd.DataFrame:
    """
    Convert a TSV or CSV file into a pandas DataFrame.
    
    Args:
        input_txt (str): TSV or CSV file containing transcript expression levels.
    
    Returns:
        df_gene (pd.DataFrame): Pandas DataFrame with 'Transcript' and 'Expression_level' as columns.
    
    Raises:
        None
    """
    df_input = pd.read_csv(
        input_txt,
        sep=r"[\t,]",
        lineterminator="\n",
        names=["Transcript", "Expression_level"],
        engine="python",
    )
    return df_input


def exprLevel_byGene(
    df_exprTranscript: pd.DataFrame, df_output_gtf_selection: pd.DataFrame
) -> pd.DataFrame:
    """
    Find the gene of each transcript given by the expression level CSV/TSV file
    and sum the expression level of all transcripts from the same gene.
    
    Args:
        df_exprTranscript (pd.DataFrame): Pandas DataFrame containing transcripts and their expression levels,
                                          generated by the "tsv_or_csv_to_df" function.
        df_output_gtf_selection (pd.DataFrame): Pandas DataFrame containing genes and transcripts,
                                                generated by the "transcripts_by_gene_inDf" function.
    
    Returns:
        Pandas DataFrame having 'Gene' and sum of its transcript expression levels.
    
    Raises:
        None
    """
    df_merged = pd.merge(df_output_gtf_selection, df_exprTranscript, how="inner", on="Transcript")
    df_sum = df_merged.groupby("Gene")["Expression_level"].sum().reset_index()
    return df_sum


def match_byGene(
    df_reprTranscript: pd.DataFrame, df_expressionLevel_byGene: pd.DataFrame
) -> pd.DataFrame:
    """
    Find matching genes between the two DataFrames.
    
    Args:
        df_reprTranscript (pd.DataFrame): Pandas DataFrame containing genes and their representative transcripts,
                                          generated by the "dict_reprTrans_to_df()" function.
        df_expressionLevel_byGene (pd.DataFrame): Pandas DataFrame containing genes and their expression levels,
                                                  generated by the "transcript_by_gene_inDf()" function.
    
    Returns:
        Pandas DataFrame having representative transcripts and their expression levels.
    
    Raises:
        None
    """
    df_merged = pd.merge(df_reprTranscript, df_expressionLevel_byGene, how="inner", on="Gene")
    df_clean = df_merged.loc[:, ["reprTrans", "Expression_level"]]
    return df_clean


# functions to run this part of the program


def match_reprTranscript_expressionLevel(
    exprTrans: str, dict_reprTrans: dict, gtf_file: str,
):
    """
    Combine functions to replace transcripts from an expression level CSV/TSV file with representative transcripts.

    Args:
        exprTrans (str): CSV or TSV file containing transcripts and their expression level.
        dict_reprTrans (dict): Dictionary of genes and their representative transcripts.
        gtf_file (str): Path to the GTF file.

    Returns:
        Pandas DataFrame of representative transcripts and their expression level.

    Raises:
        None
    """
    df_gene_transcript = gtf_to_df(gtf_file)
    df_exprTrans = tsv_or_csv_to_df(exprTrans)
    df_reprTrans = dict_reprTrans_to_df(dict_reprTrans)
    df_exprLevel_byGene = exprLevel_byGene(df_exprTrans, df_gene_transcript)
    df_match = match_byGene(df_reprTrans, df_exprLevel_byGene)
    df_match.rename(columns={"reprTrans": "id", "Expression_level": "level"}, inplace=True)
    return df_match


def transcript_sampling(total_transcript_number, df_repr, output_csv):
    total = df_repr["level"].sum()
    total_transcript_number = int(total_transcript_number)
    normalized = total_transcript_number / total
    levels = np.random.poisson(df_repr["level"] * normalized)
    transcript_numbers = pd.DataFrame({"id": df_repr["id"], "count": levels})
    transcript_numbers.to_csv(output_csv, index=False)
