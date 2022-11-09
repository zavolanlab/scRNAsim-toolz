
import pandas as pd
import re 
import itertools 

'''
This code take as input a gtf file and returns a dictionary of transcripts with best support level of each gene of the input

'''



##import modified gtf file and create a df##

def import_gtfSelection_to_df(gtf_modified_file):

    #create a df from the tab separated file input
    df_input =pd.read_csv(gtf_modified_file, sep='\t', lineterminator='\n', 
names =["Gene_mixed", "Transcript", "Support_level", "Na1", "Na2"] )

    df_input["Support_level"] = df_input["Support_level"].replace(" ", "")

    #Create a new column with only gene name from Gene_mixed column
    df_input["Gene"] = df_input["Gene_mixed"].str.extract('([A-Z]\w{0,})', expand=True)

    #Create a new column with only transcript number from Gene_mixed column
    df_input["Transcript_number"] = df_input["Gene_mixed"].str.extract('(^\d)', expand=True)

    #Create a new df with relevant column and without NA
    df_clean = df_input.loc[:, ["Gene", "Transcript","Support_level"]]
    df_clean["Gene"] = df_clean["Gene"].fillna(method='ffill')
    df_clean = df_clean.dropna(axis=0)
    return df_clean



##Returns a df containing representative transcripts and their expression level from genes mentioned in the csv file##

def representative_transcripts_inDict(df_gtfSelection): 
   

    #create a df indexed on booth Gene and Transcript columns 
    df_multIndex = df_gtfSelection.set_index(["Gene", "Transcript"])
    #create a df with only the transcripts with the highest support level (best is = 1 )
    df_min = df_multIndex.groupby(level=["Gene"])["Support_level"].transform("min")
    print("\n=== This is your 10 first representative transcripts : === \n \n {}".format(df_min.head(10)))
    #create a df without transcript levels
    df_final = df_multIndex.reset_index(level="Transcript")
    df_final = df_final.drop(columns=["Support_level"])
    
    #create a dict with only Gene and representative transcripts
    dict_representative_transcripts = df_final.groupby("Gene")["Transcript"].apply(list).to_dict()
    return dict_representative_transcripts  



### add your inputs here ! ###

gtf_file = "Homo_sapiens.GRCh38.107_intermediat_file.txt" # add the gtf input file here 

df_gtf = import_gtfSelection_to_df(gtf_file)

dictionary_of_representative_transcripts = representative_transcripts_inDict(df_gtf)
