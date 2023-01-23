import pandas as pd
import numpy as np
from gtfparse import read_gtf


def attributs_converter(attributs):
    """
    This funtion converts the "unstrucktured" ;-seperated part of he line into a list of identifyers and coresponding data the struckture of
    which can be used ot find the data easyly e.g the index of the identifier transcrip_id + 1 will give the trasncript id of the current gene
    Input: 
        attributs = str() #the unstrucktured part of the entry
    Output:
        attributs = list() # cleand list with the characterritsics discribed above
    """
    attributs = attributs.replace('"', "")
    attributs = attributs.replace(";", "")
    attributs = attributs.replace("\\n", "")
    attributs = attributs.split(" ")

    return attributs


def find_in_attributs(attributs, look_for):
    """
    This function finds a key word and used that to lokat the value of that key word e.g key = gene_id, value = 'ENSMUSG00002074970',
    this works as they are next to each other in the attributs list. 
    Inputs:
        sub_enty = list() 
        look_fore = str() #string of with the name of the key to look for
    Output: 
        attributs[index] or NA = str() #NA is returned if the key was not found in the attributs
    """
    try:
        index = attributs.index(look_for) + 1
        return attributs[index]
    except:
        # print("No",look_for,"in the entry the return was set to NA\n",attributs)
        return "NA"


def _re_format(rep_trans_dict):
    """
    This function is ment to reformat dictionary of the representatice transcripts into an dictionary with only one entry per key
    Input:
        rep_trans_dict = {gene_id : [transcript_id , transcript_support_level , transcript_length]}
    Output: 
        rep_transcripts = {gene_id : transcript_id}
    """
    rep_transcripts = dict()
    for gene_id in rep_trans_dict:
        rep_transcripts[gene_id] = rep_trans_dict[gene_id][0]

    return rep_transcripts


def get_rep_trans(file_name="test"):
    """ 
    This is the main function of this script it selects one representative transcrip per gene based on a gtf annotation file. 
    It does so be two criteria: first the transcript support level and it there are several transcript 
    of one gene that have the same trasncript_support_level it chooses the one that corresponds to the longest mRNA.
    Input: 
        file_name = str() # name of the annotation file with or without the .gtf part
    Output: 
        rep_transcripts = {gene_id : transcript_id}
    """

    # setting defoult variables
    rep_trans = dict()
    cur_gID = str()
    cur_best_trans = [
        str(),
        100,
        0,
    ]  # [transcript_id , transcript_support_level , transcript_length]
    pot_best_trans = False
    cur_tID = str()
    ignor_trans = False

    with open(file_name, "r") as f:
        for line in f:
            entry = line.split("\t")

            # removes expected but unneeded entrys
            exp_unneed = [
                "CDS",
                "stop_codon",
                "five_prime_utr",
                "three_prime_utr",
                "start_codon",
                "Selenocysteine",
            ]
            if len(entry) == 1 or entry[2] in exp_unneed:
                continue

            # this function turns the less organized part of the entry into a reable list
            attributs = attributs_converter(entry[8])
            # looking for and processing exons entrys
            if entry[2] == "exon":

                # dicide if to contiune or not
                if ignor_trans:
                    continue
                elif cur_gID != attributs[1]:
                    raise ValueError("ERROR exon from an unexpected Gen")
                    continue
                elif find_in_attributs(attributs, "transcript_id") != cur_tID:
                    raise ValueError("exon from an unexpected transcript")
                    continue

                # adding the length of the exon to the appropriat list and chacking for changes in best transcript
                if pot_best_trans:
                    pot_best_trans[2] += int(entry[4]) - int(entry[3])
                    if pot_best_trans[2] > cur_best_trans[2]:
                        cur_best_trans = pot_best_trans
                        pot_best_trans = False
                else:
                    cur_best_trans[2] += int(entry[4]) - int(entry[3])

            # looking for and processing transcript entrys
            elif entry[2] == "transcript":

                # varryfi that the gen is correct
                if cur_gID != attributs[1]:
                    raise ValueError("ERROR transcript from an unexpected Gen")
                    continue

                # finding the transcript id and the support level
                cur_tID = find_in_attributs(attributs, "transcript_id")
                t_supp_lvl = find_in_attributs(attributs, "transcript_support_level")

                # If there is no transcript support level or the level is given as NA it is nomed as 100. else the transcript support level is tunrn into int
                if t_supp_lvl == "NA":
                    t_supp_lvl = 100
                else:
                    try:
                        t_supp_lvl = int(t_supp_lvl)
                    except:
                        t_supp_lvl = 100

                # decides if the transcript has potential to become the representative transcript
                if t_supp_lvl < cur_best_trans[1] or cur_best_trans[0] == "":
                    cur_best_trans = [cur_tID, t_supp_lvl, 0]
                    pot_best_trans = False
                    ignor_trans = False

                elif t_supp_lvl == cur_best_trans[1]:
                    pot_best_trans = [cur_tID, t_supp_lvl, 0]
                else:
                    ignor_trans = True

            # looking for and processing gene entrys
            elif entry[2] == "gene":

                # updating rep_trans dict
                if cur_gID not in rep_trans:
                    rep_trans[cur_gID] = cur_best_trans
                else:
                    if rep_trans[cur_gID][1] > cur_best_trans[1]:
                        rep_trans[cur_gID] = cur_best_trans
                    elif (
                        rep_trans[cur_gID][1] == cur_best_trans[1]
                        and rep_trans[cur_gID][2] < cur_best_trans[2]
                    ):
                        rep_trans[cur_gID] = cur_best_trans

                # updating cur_gID and resetting cur_best_trans
                cur_gID = attributs[1]
                cur_best_trans = [str(), 100, 0]

            # raises an error for unidentifyable entrys
            else:
                raise ValueError("This entry could not be identified\n", entry)

        # addding the final gene to the dictionary
        if cur_gID not in rep_trans:
            rep_trans[cur_gID] = cur_best_trans
        else:
            if rep_trans[cur_gID][1] > cur_best_trans[1]:
                rep_trans[cur_gID] = cur_best_trans
            elif (
                rep_trans[cur_gID][1] == cur_best_trans[1]
                and rep_trans[cur_gID][2] < cur_best_trans[2]
            ):
                rep_trans[cur_gID] = cur_best_trans

        del rep_trans[""]
        rep_transcripts = _re_format(rep_trans)
        return rep_transcripts


def _test():
    """
    This funtion is ment to be run for test
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
        print("The test was succses full")


def gtf_file_writer(original_file, rep_transcript_dict, output_file):
    """
    this function writes the output GTF file
    """
    output = []

    with open(original_file, "r") as f:
        for line in f:
            entry = line.split("\t")
            if line[0] != "#":
                attributes = attributs_converter(entry[8])
                type_ = entry[2]
            else:
                continue
            if type_ == "gene":
                gene_id = find_in_attributs(attributes, "gene_id")
                output.append(line)
            else:
                transcript_id = find_in_attributs(attributes, "transcript_id")
                if rep_transcript_dict[gene_id] == transcript_id:
                    output.append(line)

    with open(output_file, "w") as last_file:
        for item in output:
            last_file.write(item)


def gtf_to_df(gtf_file: str) -> pd.DataFrame:
    """
    This function take a .gtf file and convert it into a 
    dataframe containing gene_id and their transcripts_id.
        Args:
            gtf_file (str) : path to the .gtf file

        Returns:
            df_gtf (pd.DataFrame) : pandas dataframe containing columns
            gene_id and their transcripts_id.
        Raises : 
            None 
    
    """
    df_gtf = read_gtf(gtf_file)
    df_gtf = df_gtf.loc[df_gtf["feature"] == "transcript"]
    df_gtf = df_gtf[["gene_id", "transcript_id"]]
    df_gtf = df_gtf.rename(columns={"gene_id": "Gene", "transcript_id": "Transcript"})
    return df_gtf


def dict_reprTrans_to_df(dict_reprTrans: dict[str, str]) -> pd.DataFrame:

    """Convert a dictionary of genes and their representative transcript into a dataframe 

        Args:
            dict_reprTrans (dict) : {'Gene':['transcriptA', 'transcriptB'], ...}

        Returns:
            Pandas dataframe having Gene and transcript as columns
      
        Raises:
            Only dict are allowed
            Key should be strings
            Value should be strings
          
    """
    pass
    if not type(dict_reprTrans) is dict:
        raise TypeError("Only dict are allowed")
    if type(list(dict_reprTrans.keys())[0]) is not str:
        raise TypeError("Key should be strings")
    if type(list(dict_reprTrans.values())[0]) is not str:
        raise TypeError("Values should be strings")

    df_reprTrans = pd.DataFrame.from_dict(
        dict_reprTrans, orient="index", columns=["reprTranscript"]
    )
    df_reprTrans = df_reprTrans.reset_index(level=0)
    df_reprTrans.columns = ["Gene", "reprTrans"]
    df_reprTrans["reprTrans"] = df_reprTrans["reprTrans"].str.replace(
        r"\.[1-9]", "", regex=True
    )
    return df_reprTrans


def tsv_or_csv_to_df(input_txt: str) -> pd.DataFrame:
    """Convert tsv or csv file into a pandas dataframe

        Args:
            input_txt (str): csv or tsv file containing transcript expression level

        Returns:
            df_gene (str): Pandas dataframe having transcript and expression level
            as columns  
      
        Raises:
            None          
    """
    pass
    df_input = pd.read_csv(
        input_txt,
        sep=r"[\t,]",
        lineterminator="\n",
        names=["Transcript", "Expression_level"],
        engine="python",
    )
    return df_input


def exprLevel_byGene(
    df_exprTrasncript: pd.DataFrame, df_output_gtf_selection: pd.DataFrame
) -> pd.DataFrame:
    """find the gene of each transcipt given by the expression level csv/tsv file,
       and summ expression level of all transcipts from the same gene. 

        Args:
            df_exprTranscript : pandas Dataframe containing transcript and their expression level,
            generated by "tsv_or_csv_to_df" function
            df_output_gtf_selection : pandas Dataframe containing genes and transcripts,
            generated by "transcripts_by_gene_inDf" function 

        Returns:
            Pandas dataframe having gene and sum of its transcript expression level
      
        Raises:
            None          
    """
    pass
    df_merged = pd.merge(
        df_output_gtf_selection, df_exprTrasncript, how="inner", on="Transcript"
    )
    df_sum = df_merged.groupby("Gene").sum("Expression_level")
    return df_sum


def match_byGene(
    df_reprTranscript: pd.DataFrame, df_expressionLevel_byGene: pd.DataFrame
) -> pd.DataFrame:
    """Find matching genes bewteen the 2 args 

        Args:
            df_reprTranscript : pandas Dataframe containing genes 
            and their representative transcript, generated by
            "dict_reprTrans_to_df()" 
            df_expressionLevel_byGene : pandas Dataframe containing 
            genes and their expression level generated by 
            "transcript_by_gene_inDf()"

        Returns:
            Pandas dataframe having representative trasncripts 
            and their expression level
      
        Raises:
            None          
    """
    pass
    df_merged = pd.merge(
        df_reprTranscript, df_expressionLevel_byGene, how="outer", on="Gene"
    )
    df_clean = df_merged.dropna(axis=0)
    df_clean = df_clean.loc[:, ["reprTrans", "Expression_level"]]
    return df_clean


### functions to run this part of the programm


def match_reprTranscript_expressionLevel(
    exprTrans: str, dict_reprTrans: dict, gtf_file: str,
):
    """Combine functions to replace transcripts from an expression level csv/tsv file 
       with representative transcripts 

        Args:
            exprTrans (str): csv or tsv file containing transcripts
            and their expression level 
            dict_reprTrans (dict) : dict of genes and their 
            representative transcipt
            intemediate_file (str) : txt file containing genes, transcript 
            and their expression level from the transkript_extractor function
            output_path : path indicating were the tsv file should be written

        Returns:
            tsv file of representative trasncripts and their expression level
      
        Raises:
            None          
    """
    df_gene_transcript = gtf_to_df(gtf_file)
    df_exprTrans = tsv_or_csv_to_df(exprTrans)
    df_reprTrans = dict_reprTrans_to_df(dict_reprTrans)
    df_exprLevel_byGene = exprLevel_byGene(df_exprTrans, df_gene_transcript)
    df_match = match_byGene(df_reprTrans, df_exprLevel_byGene)
    df_match.rename(
        columns={"reprTrans": "id", "Expression_level": "level"}, inplace=True
    )
    return df_match


def transcript_sampling(total_transcript_number, df_repr, output_csv):
    df = df_repr
    levels = []
    sums = df["level"].tolist()
    total = sum(sums)
    total_transcript_number = int(total_transcript_number)
    normalized = total_transcript_number / total
    for expression_level in df["level"]:
        poisson_sampled = np.random.poisson(expression_level * normalized)
        levels.append(poisson_sampled)

    transcript_numbers = pd.DataFrame({"id": df["id"], "count": levels})
    pd.DataFrame.to_csv(transcript_numbers, output_csv)
