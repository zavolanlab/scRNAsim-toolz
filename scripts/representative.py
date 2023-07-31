'''
This part of the code take as input a gtf modified file
and return a dictionary of transcripts with best
support level for each gene of the input
'''
import pandas as pd
# import os


def import_gtf_selection_to_df(gtf_modified_file: str) -> pd.DataFrame:
    """Import intermediate file from gtf and create a df

        Args:
            gtf_modified_file (str) : path to the intermediate file

        Returns:
            Pandas dataframe having Gene, transcript
            and support level as columns

        Raises:
            TypeError : Only str path is allowed

    """
    if not isinstance(gtf_modified_file, str):
        raise TypeError("Only str path is allowed")
    df_input = pd.read_csv(
        gtf_modified_file, sep='\t', lineterminator='\n',
        names=["Gene_mixed", "Transcript", "Support_level", "Na1", "Na2"]
        )
    df_input["Support_level"] = df_input["Support_level"].replace(" ", "")
    df_input["Gene"] = df_input["Gene_mixed"].str.extract(
        r'([A-Z]\w{0,})', expand=True  # noqa: W605
        )
    df_input["Transcript_number"] = df_input["Gene_mixed"].str.extract(
        r'(^\d)', expand=True  # noqa: W605
        )
    df_clean = df_input.loc[:, ["Gene", "Transcript", "Support_level"]]
    df_clean["Gene"] = df_clean["Gene"].fillna(method='ffill')
    df_clean = df_clean.dropna(axis=0)
    return df_clean


def representative_transcripts_in_dict(
        df_gtf_selection: pd.DataFrame) -> pd.DataFrame:
    """Return a dict containing for each gene transcripts
        with highest confidence level

        Args:
            df_gtf_selection (str): Pandas dataframe having Gene,
            transcript and support level as columns

        Returns:
            Dict {'Gene':['transcriptA', 'transcriptB'], ...}

        Raises:
            TypeError : Only pandas DataFrame is allowed
    """
    if not isinstance(df_gtf_selection, pd.DataFrame):
        raise TypeError("Only pandas DataFrame is allowed")
    df_min = df_gtf_selection[
        df_gtf_selection["Support_level"] ==
        df_gtf_selection.groupby("Gene")["Support_level"].transform(min)
        ]
    df_final = df_min.drop(columns=["Support_level"])
    dict_representative_transcripts = df_final.groupby("Gene")[
        "Transcript"].apply(list).to_dict()
    return dict_representative_transcripts


def find_repr_by_support_level(intermediate_file: str) -> dict[str, str]:
    """Combine functions import_gtf_selection_to_df()
        and representative_transcripts_in_dict()

        Args:
            intermediate_file : path to the intermediate file

        Returns:
            Dict {'Gene':['transcriptA', 'transcriptB'], ...}

        Raises:
            None


    """
    df_gtf = import_gtf_selection_to_df(intermediate_file)
    dict_repr_trans = representative_transcripts_in_dict(df_gtf)
    return dict_repr_trans


# if __name__ == "__main__":
#     find_repr_by_support_level()
