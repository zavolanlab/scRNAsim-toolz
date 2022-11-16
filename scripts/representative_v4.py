
import pandas as pd

'''
This part of the code take as input a gtf modified file 
and return a dictionary of transcripts with best
support level for each gene of the input

'''




def import_gtfSelection_to_df(gtf_modified_file: str):
    """Import intermediate file from gtf and create a df

        Args:
            gtf_modified_file (str) : path to the intermediate file

        Returns:
            Pandas dataframe having Gene, transcript 
            and support level as columns
      
        Raises:
            TypeError : Only str path is allowed
          
    """
    pass
    if not type(gtf_modified_file) is str:
      raise TypeError("Only str path is allowed")
    df_input = pd.read_csv(gtf_modified_file, sep = '\t', lineterminator = '\n', 
names = ["Gene_mixed", "Transcript", "Support_level", "Na1", "Na2"] )
    df_input["Support_level"] = df_input["Support_level"].replace(" ", "")
    df_input["Gene"] = df_input["Gene_mixed"].str.extract('([A-Z]\w{0,})', expand=True)
    df_input["Transcript_number"] = df_input["Gene_mixed"].str.extract('(^\d)', expand=True)
    df_clean = df_input.loc[:, ["Gene", "Transcript","Support_level"]]
    df_clean["Gene"] = df_clean["Gene"].fillna(method = 'ffill')
    df_clean = df_clean.dropna(axis = 0)
    return df_clean




def representative_transcripts_inDict(df_gtfSelection: str) -> pd.DataFrame:
    """Return a dict containing for each gene transcripts 
        with highest confidence level

        Args:
            df_gtfSelection (str): Pandas dataframe having Gene,
            transcript and support level as columns

        Returns:
            Dict {'Gene':['transcriptA', 'transcriptB'], ...}
      
        Raises:
            TypeError : Only pandas DataFrame is allowed
    """
    pass 

    if not type(df_gtfSelection) is pd.DataFrame:
        raise TypeError("Only pandas DataFrame is allowed")
  
    df_multIndex = df_gtfSelection.set_index(["Gene", "Transcript"])
    #highest support level = 1 , worst = 5, NA = 100
    df_min = df_multIndex.groupby(level=["Gene"])["Support_level"].transform("min")
    df_final = df_min.reset_index(level = "Transcript")
    df_final = df_final.drop(columns = ["Support_level"])
    dict_representative_transcripts = df_final.groupby("Gene")["Transcript"].apply(list).to_dict()
    return dict_representative_transcripts  



def find_repr_by_SupportLevel(intermediate_file:str): 
    """Combine functions import_gtfSelection_to_df() 
        and representative_transcripts_inDict()

        Args:
            intermediate_file : path to the intermediate file

        Returns:
            Dict {'Gene':['transcriptA', 'transcriptB'], ...}
      
        Raises:
            None

          
    """
    pass 
    df_gtf = import_gtfSelection_to_df(intermediate_file)
    dict_reprTrans = representative_transcripts_inDict(df_gtf)
    return dict_reprTrans



if __name__ == "__main__":  
    find_repr_by_SupportLevel() 
