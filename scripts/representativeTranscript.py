
import pandas as pd
import re 

'''
This part of the code find a representative transcript  for each ID given in the csv input file. 
It works with the output of the function "transcript selector".
It return a csv file containing the names of the representative transcript and their levels.

'''

#import csv file [ID, Expression_level] and create a df
def import_csv_to_df(csv_file) :
    df_csv = pd.read_csv(csv_file, names=["ID", "Expression_level"])
    return df_csv


#import modified gtf file and create a df
def import_gtfSelection_to_df(gtf_modified_file):

    #create a df from the tab separated file input
    df_input =pd.read_csv(gtf_modified_file, sep='\t', lineterminator='\n', 
names =["Gene_mixed", "Transcript", "Support_level", "Length", "Line_index", "Na"] )

    #Create a new column with only gene name from Gene_mixed column
    df_input["Gene"] = df_input["Gene_mixed"].str.extract('([A-Z]\w{0,})', expand=True)

    #Create a new column with only transcript number from Gene_mixed column
    df_input["Transcript_number"] = df_input["Gene_mixed"].str.extract('(^\d)', expand=True)

    #Create a new df with relevant column and without NA
    df_clean = df_input.loc[:, ["Gene","Transcript_number","Transcript","Support_level", "Length", "Line_index"]]
    df_clean["Gene"] = df_clean["Gene"].fillna(method='ffill')
    df_clean = df_clean.dropna(axis=0)
    return df_clean



# Return a df containing representative transcripts and their expression level from genes mentioned in the csv file 

def representative_from_gene(df_gtfSelection, df_csv): 
    #check wich genes are shared by df_gtfSelection and df_csv and create a reduced df
    df_shared = df_gtfSelection[(df_gtfSelection["Gene"].isin(df_csv["ID"]))]# /!\ return an empty list if I do the same with "Transcript" column. 

    #pick only the transcripts with the highest support level (best is = 1 )
    df_filtered = df_shared[df_shared["Support_level"]==df_shared["Support_level"].min()]

    #pick the transcript with the greatest length if multiple transcript have the same support level
    df_filtered2 = df_filtered.sort_index().groupby(df_filtered["Gene"]).max()

    #combine expression level with transcript with highest support level and greatest length 
    df_filtered2["Expression_level"] = df_filtered2.Gene.map(df_csv.set_index("ID")['Expression_level'])

    #create a reduced df with only the representative transcript and its support level 
    df_final = df_filtered2[["Transcript","Expression_level"]]
    print("This is your representative transcripts : \n \n {}".format(df_final))

    return df_final 


#Missing conditions :
    # If multiple transcript form the same gene, add their expression 
    # If given transcript name, return the most representative transcript 



# Output csv file containing representative transcript levels and id
def write_csv_Transcript_Expression_level(df_representative):
    df_final = df_representative[["Transcript","Expression_level"]]
    with open("transcripts.csv", "w") as fileout : 
        fileout.write(df_representative.to_csv(columns=["Transcript","Expression_level"], index=False))

### put your inputs here ! ###

ID_levels_csv = "scRNA\input.csv" # put your csv file here
gtf_file = "scRNA/test_inter_mediat_file.txt" # put your gtf selected file here 

df_gtf = import_gtfSelection_to_df(gtf_file)
df_csv_input = import_csv_to_df(ID_levels_csv)

representative_transcript = representative_from_gene(df_gtf, df_csv_input)

write_csv_Transcript_Expression_level(representative_transcript) # create a csv file named "transcript.csv"

