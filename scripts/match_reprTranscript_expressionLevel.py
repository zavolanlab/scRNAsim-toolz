
import pandas as pd
import json
import re
import rerpresentative_v4 as repr
import os


def dict_reprTrans_to_df(dict_reprTrans: dict):

    """Convert a dictionary of genes and their representative transcript into a dataframe 

        Args:
            dict_reprTrans (dict) : {'Gene':['transcriptA', 'transcriptB'], ...}

        Returns:
            Pandas dataframe having Gene and transcript as columns
      
        Raises:
            /!\ None, I wasn't able to make a TypeError with dict  
            : Only dict made of key string and value string is allowed
          
    """
    pass

    df_reprTrans = pd.DataFrame.from_dict(dict_reprTrans, orient="index", columns=["reprTranscript"])
    df_reprTrans = df_reprTrans.reset_index(level=0)
    df_reprTrans.columns = ["Gene", 'reprTrans']
    df_reprTrans["reprTrans"] = df_reprTrans["reprTrans"].str.replace(r'\.[1-9]', '', regex=True)
    return df_reprTrans


def txt_to_dict(dict_txt: str):
    """Convert a txt file into a dictionary 

        Args:
            dict_txt (str) : pathe to a txt file of a dict
            structured as {'Gene':['transcriptA', 'transcriptB'], ...}

        Returns:
            dict (dict) : dictionary stuctured as {'Gene':['transcriptA', 'transcriptB'], ...}
      
        Raises:
            None          
    """
    pass

    input : str = open(dict_txt, "r").read()
    input : str = input.replace("\'", "\"")
    dict = json.loads(input)
    return dict



def transcripts_by_gene_inDf(df_gtfSelection: str) -> pd.DataFrame:
    """Convert multiindex dataframe from function into a simple dataframe 

        Args:
            df_gtfSelection (str): Pandas multiindex dataframe having Gene,
            transcript as indexs and support level as columns. 
            Come from the function import_gtfSelection_to_df()

        Returns:
            df_gene (str): Pandas dataframe having Gene and
            transcript as columns 
      
        Raises:
            None          
    """
    pass
    df_gene = df_gtfSelection.set_index(["Gene"])
    df_gene = df_gene.drop(columns=["Support_level"])
    df_gene['Transcript']=df_gene['Transcript'].str.replace(r"\.[0-9]","", regex=True)
    df_gene = df_gene.reset_index(level=0)
    return df_gene


def tsv_or_csv_to_df(input_txt:str) :
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
    df_input =pd.read_csv(input_txt, sep=r"[\t,]", lineterminator='\n',
     names=["Transcript", "Expression_level"],
     engine = "python")
    return df_input


def exprLevel_byGene(df_exprTrasncript:str, df_output_gtf_selection:str) -> pd.DataFrame :
    """Find matching transcripts bewteen the 2 args 

        Args:
            df_exprTranscript (str): pandas Dataframe containing transcript and their expression level
            df_output_gtf_selection (str) : pandas Dataframe containing genes and transcripts 

        Returns:
            Pandas dataframe having gene and sum of its transcript expression level
      
        Raises:
            None          
    """
    pass 
    df_merged = pd.merge(df_output_gtf_selection, df_exprTrasncript , how="inner", on="Transcript")
    df_sum = df_merged.groupby("Gene").sum("Expression_level") # sum transcripts comming from the same gene  
    return df_sum

def match_byGene(df_reprTranscript:str, df_expressionLevel_byGene:str) -> pd.DataFrame: 
    """Find matching genes bewteen the 2 args 

        Args:
            df_reprTranscript (str): pandas Dataframe containing genes 
            and their representative transcript
            df_expressionLevel_byGene (str) : pandas Dataframe containing 
            genes and their expression level 

        Returns:
            Pandas dataframe having representative trasncripts 
            and their expression level
      
        Raises:
            None          
    """
    pass 
    df_merged = pd.merge(df_reprTranscript, df_expressionLevel_byGene , how="outer", on="Gene")
    df_clean = df_merged.dropna(axis=0)
    df_clean = df_clean.loc[:, ["reprTrans","Expression_level"]]
    return df_clean

def output_tsv(dataframe:str)-> pd.DataFrame :
    """Convert pandas dataframe into a tsv file 

        Args:
            dataframe (str): Pandas dataframe containing
            representative transcripts and their expression level 

        Returns:
            Tsv file containing representative transcripts
             and their expression level in the same directory
      
        Raises:
            None          
    """
    pass 

    csv_file = dataframe.to_csv(os.getcwd()+"\ReprTrans_ExpressionLevel.tsv", sep="\t", 
    index=False, header=False)
    return csv_file

### functions to run this part of the programm

def match_reprTranscript_expressionLevel(exprTrans:str, dict_reprTrans:dict, intermediate_file:str): 
    """Combine functions to replace transcripts from an expression level csv/tsv file 
       with representative transcripts 

        Args:
            exprTrans (str): csv or tsv file containing transcripts
            and their expression level 
            dict_reprTrans (dict) : dict of genes and their 
            representative transcipt
            intemediate_file (str) : txt file containing genes, transcript 
            and their expression level from the transkript_extractor function

        Returns:
            tsv file of representative trasncripts and their expression level
      
        Raises:
            None          
    """
    df_intermediate = repr.import_gtfSelection_to_df(intermediate_file)
    df_geneTrans = transcripts_by_gene_inDf(df_intermediate)
    df_exprTrans = tsv_or_csv_to_df(exprTrans)
    df_reprTrans = dict_reprTrans_to_df(dict_reprTrans)
    df_exprLevel_byGene = exprLevel_byGene(df_exprTrans, df_geneTrans)
    df_match = match_byGene(df_reprTrans, df_exprLevel_byGene)
    output = output_tsv(df_match)
    return output


# run the programm 

dict_txt = a #input a dict of {gene:reprTrans} in the form of a txt file
input_intermediate_file = b #input the intermediate file generated by transckript extractor
input_expr = c #input a csv or tsv file containing the expr level 

dict_reprTrans = txt_to_dict(dict_txt)
match_final = match_reprTranscript_expressionLevel(input_expr, dict_reprTrans, input_intermediate_file)
print("this is the function :\n\n {}".format(match_final))

if __name__ == "__main__":  
    match_reprTranscript_expressionLevel()
 