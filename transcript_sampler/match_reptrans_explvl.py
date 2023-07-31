"""Match representative transcript with expression level"""
# Made by Hugo Gillet #

import logging
import pandas as pd
from gtfparse import read_gtf

LOG = logging.getLogger(__name__)


class MatchReptransExplvl:
    """Match representative transcript with expression level"""
    def __init__(self):
        pass

    @staticmethod
    def gtf_to_df(gtf_file: str) -> pd.DataFrame:
        """
        This function takes a .gtf file and converts it into a pandas DataFrame
        containing gene_id and their transcript_id.

        Args:
            gtf_file (str): Path to the .gtf file.

        Returns:
            df_gtf (pd.DataFrame): Pandas DataFrame containing columns
            'Gene' and 'Transcript'.

        Raises:
            None
        """
        df_gtf = read_gtf(gtf_file,).to_pandas()
        df_gtf = df_gtf[df_gtf["feature"] == "transcript"]
        df_gtf = df_gtf[["gene_id", "transcript_id"]]
        df_gtf = df_gtf.rename(columns={
            "gene_id": "Gene", "transcript_id": "Transcript"
            })
        return df_gtf

    @staticmethod
    def dict_repr_trans_to_df(dict_reprTrans: "dict[str, str]") -> pd.DataFrame:
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
            LOG.error("Only dictionaries are allowed")
            raise TypeError("Only dictionaries are allowed")
        if not all(isinstance(key, str) for key in dict_reprTrans.keys()):
            LOG.error("Keys should be strings")
            raise TypeError("Keys should be strings")
        if not all(isinstance(value, str) for value in dict_reprTrans.values()):
            LOG.error("Values should be strings")
            raise TypeError("Values should be strings")

        df_reprTrans = pd.DataFrame.from_dict(dict_reprTrans, orient="index", columns=["reprTranscript"])
        df_reprTrans = df_reprTrans.reset_index()
        df_reprTrans.columns = ["Gene", "reprTrans"]
        df_reprTrans["reprTrans"] = df_reprTrans["reprTrans"].str.replace(r"\.[1-9]", "", regex=True)

        return df_reprTrans

    @staticmethod
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

    @staticmethod
    def expr_level_by_gene(
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

    @staticmethod
    def match_by_gene(
        df_reprTranscript: pd.DataFrame, df_expressionLevel_byGene: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Find matching genes between the two DataFrames.
        
        Args:
            df_reprTranscript (pd.DataFrame): Pandas DataFrame containing genes and their representative transcripts,
                                            generated by the "dict_repr_trans_to_df()" function.
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
    
    def match_repr_transcript_expression_level(
        self, exprTrans: str, dict_reprTrans: dict, gtf_file: str,
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
        df_gene_transcript = self.gtf_to_df(gtf_file)
        df_exprTrans = self.tsv_or_csv_to_df(exprTrans)
        df_reprTrans = self.dict_repr_trans_to_df(dict_reprTrans)
        df_expr_level_by_gene = self.expr_level_by_gene(df_exprTrans, df_gene_transcript)
        df_match = self.match_by_gene(df_reprTrans, df_expr_level_by_gene)
        df_match.rename(columns={"reprTrans": "id", "Expression_level": "level"}, inplace=True)
        return df_match



# def dict_repr_trans_to_df(dict_reprTrans: "dict[str, str]") -> pd.DataFrame:

#     """Convert a dictionary of genes and their representative
#     transcript into a dataframe

#         Args:
#             dict_reprTrans (dict): {'Gene':['transcriptA', 'transcriptB'], ...}

#         Returns:
#             Pandas dataframe having Gene and transcript as columns

#         Raises:
#             Only dict are allowed
#             Key should be strings
#             Value should be strings

#     """
#     pass
#     if not type(dict_reprTrans) is dict:
#         raise TypeError("Only dict are allowed")
#     if type(list(dict_reprTrans.keys())[0]) is not str:
#         raise TypeError("Key should be strings")
#     if type(list(dict_reprTrans.values())[0]) is not str:
#         raise TypeError("Values should be strings")

#     df_reprTrans = pd.DataFrame.from_dict(
#         dict_reprTrans, orient="index", columns=["reprTranscript"]
#     )
#     df_reprTrans = df_reprTrans.reset_index(level=0)
#     df_reprTrans.columns = ["Gene", "reprTrans"]
#     df_reprTrans["reprTrans"] = df_reprTrans["reprTrans"].str.replace(
#         r"\.[1-9]", "", regex=True
#     )
#     return df_reprTrans


# def gene_and_transcript(gtf_file: str) -> pd.DataFrame:
#     """
#     This function take a .gtf file and convert it into a
#     dataframe containing gene_id and their transcripts_id.
#         Args:
#             gtf_file(str) : path to the .gtf file

#         Returns:
#             df_gtf(pd.DataFrame): pandas df containing having has columns
#             gene_id and their transcripts_id.
#         Raises:
#             None
#     """
#     df_gtf = read_gtf(gtf_file)
#     df_gtf = df_gtf.loc[df_gtf["feature"] == "transcript"]
#     df_gtf = df_gtf[["gene_id", "transcript_id"]]
#     df_gtf = df_gtf.rename(columns={"gene_id": "Gene",
#                                     "transcript_id": "Transcript"})
#     return df_gtf


# def tsv_or_csv_to_df(input_txt: str) -> pd.DataFrame:
#     """Convert tsv or csv file into a pandas dataframe

#         Args:
#             input_txt (str): csv or tsv file containing transcript exp level

#         Returns:
#             df_gene (str): Pandas dataframe having transcript and exp level
#             as columns

#         Raises:
#             None
#     """
#     pass
#     df_input = pd.read_csv(
#         input_txt,
#         sep=r"[\t,]",
#         lineterminator="\n",
#         names=["Transcript", "Expression_level"],
#         engine="python",
#     )
#     return df_input


# def expr_level_by_gene(
#     df_exprTrasncript: pd.DataFrame, df_output_gtf_selection: pd.DataFrame
# ) -> pd.DataFrame:
#     """find the gene of each transcipt given by the expression level csv/tsv
#     file, and summ expression level of all transcipts from the same gene.

#         Args:
#             df_exprTranscript: pandas df containing transcript and
#             their exp level generated by "tsv_or_csv_to_df" function
#             df_output_gtf_selection : pandas df containing genes and
#             transcripts, generated by "transcripts_by_gene_inDf" function

#         Returns:
#             Pandas dataframe having gene and sum of its transcript exp level

#         Raises:
#             None
#     """
#     pass
#     df_merged = pd.merge(
#         df_output_gtf_selection, df_exprTrasncript,
#         how="inner", on="Transcript"
#     )
#     df_sum = df_merged.groupby("Gene").sum(
#         "Expression_level"
#     )
#     return df_sum


# def match_by_gene(
#     df_reprTranscript: pd.DataFrame, df_expressionLevel_byGene: pd.DataFrame
# ) -> pd.DataFrame:
#     """Find matching genes bewteen the 2 args

#         Args:
#             df_reprTranscript : pandas Dataframe containing genes
#             and their representative transcript, generated by
#             "dict_repr_trans_to_df()"
#             df_expressionLevel_byGene : pandas Dataframe containing
#             genes and their expression level generated by
#             "transcript_by_gene_inDf()"

#         Returns:
#             Pandas dataframe having representative trasncripts
#             and their expression level

#         Raises:
#             None
#     """
#     pass
#     df_merged = pd.merge(
#         df_reprTranscript, df_expressionLevel_byGene, how="outer", on="Gene"
#     )
#     df_clean = df_merged.dropna(axis=0)
#     df_clean = df_clean.loc[:, ["reprTrans", "Expression_level"]]
#     return df_clean


# # functions to run this part of the programm
# def match_repr_transcript_expression_level(
#     exprTrans: str, dict_reprTrans: dict, gtf_file: str,
# ):
#     """Combine functions to replace transcripts from an exp level csv/tsv file
#        with representative transcripts

#         Args:
#             exprTrans (str): csv or tsv file containing transcripts
#             and their expression level
#             dict_reprTrans (dict) : dict of genes and their
#             representative transcipt
#             intemediate_file (str) : txt file containing genes, transcript
#             and their expression level from the transkript_extractor function
#             output_path : path indicating were the tsv file should be written

#         Returns:
#             tsv file of representative trasncripts and their expression level

#         Raises:
#             None
#     """
#     df_gene_transcript = gene_and_transcript(gtf_file)
#     df_exprTrans = tsv_or_csv_to_df(exprTrans)
#     df_reprTrans = dict_repr_trans_to_df(dict_reprTrans)
#     df_expr_level_by_gene = expr_level_by_gene(
#         df_exprTrans, df_gene_transcript
#         )  # error here
#     df_match = match_by_gene(df_reprTrans, df_expr_level_by_gene)
#     df_match.rename(columns={'reprTrans': 'id', 'Expression_level': 'level'},
#                     inplace=True)
#     return df_match


# # run the program
# if __name__ == "__main__":
#     match_repr_transcript_expression_level()
