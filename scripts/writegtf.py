### Called Packages ###
import pandas as pd
import numpy as np
import argparse
import re

import transcript_extractor as te

python_version = "3.7.13"
module_list =[pd,np,argparse,re]
modul_name_list = ["pd","np","argparse","re"]
### Functions ###

def transcript_ID_finder (entry):
    index_transcript_id = entry.find("transcript_id")
    find_transcript_id_name = re.compile("\"\S{1,25}\"")
    sub_entry = entry[index_transcript_id:]
    try_find_transcript_id_name = find_transcript_id_name.search(sub_entry)   
    
    try: 
        transcript_ID = try_find_transcript_id_name[0].replace("\"","")
    except:
        transcript_ID = ""
    return (transcript_ID)


'''gtf_file_writer takes as input the original gtf file and the csv file containing relevant transcripts.

    It produces a gtf file containing only the transcript entries of those contained in the csv file
    
    based on id'''


def gtf_file_writer (original_file, csv_file, output_file): 
    output = []

    df = pd.read_csv(csv_file)
    listoftranscripts = df['id'].tolist()
    if df['id'].empty:
        print('Error. \'id\' column needed in input csv file.')

    with open(original_file, 'r') as f:
            for entry in f: 
                if "\ttranscript\t" in entry:
                    transcript_id = transcript_ID_finder(entry)
                    if transcript_id in listoftranscripts:
                        output.append(entry)
    with open(output_file, 'w') as last_file:
        for line in output : # I had to add this loop because I had an error about you cannot write list in directly in a file
            last_file.write(line)


if __name__ == '__main__':
    te.version_control(module_list,modul_name_list,python_version)
    parser = argparse.ArgumentParser(
        description="gtf output file writer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--annotation", required=True, help="gtf file with genome annotation")
    parser.add_argument("--output_gtf", required=True, help="output gtf file")
    parser.add_argument("--input_csv", required=True, help="input csv file")
    args = parser.parse_args()

    gtf_file_writer(args.annotation, args.input_csv, args.output_gtf)
