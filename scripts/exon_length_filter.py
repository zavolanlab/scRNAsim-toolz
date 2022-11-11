#### Exon length filter #####

### Called Packages ###
import re
import os

import transcript_extractor as te
### Functions ###

def exon_length_calculator(entry): 
    """This funtion finds the start and end cordinates of the exon and uses them to calculate its lenght"""
    try:
        find_exon_coordinates = re.compile("\t\d{1,15}\t")
        try_find_start_coordinates = find_exon_coordinates.search(entry)
        start_coordinates = int(try_find_start_coordinates[0].replace("\t",""))
        final_index_start_coordinates = entry.find(try_find_start_coordinates[0])+len(try_find_start_coordinates[0])-1  
        sub_entry = entry[final_index_start_coordinates:]
        try_find_end_coordinates = find_exon_coordinates.search(sub_entry)
        end_coordinates = int(try_find_end_coordinates[0].replace("\t",""))
        exon_length = end_coordinates-start_coordinates
    except:
        print("\n\nIn the following enty only one or no valid coordinates could be found:\n",entry,"the value will be set to NA")
        exon_length = "NA"
    return(exon_length)

def exon_fider(entry):
    """This funtion determines if a given entry belongs to an exon"""
    exon_test = entry.find("\texon\t")
    if exon_test == -1: 
        try_exon_test = False
    else:
        try_exon_test = True
    return(try_exon_test)

def __longest_transcript_finder(current_exon_length,longest_transcript,longest_transcript_ID,old_transcript_ID):
    """This function encapsulates an operation that has to be carried out multiple times in the exon_length_filter"""
    if current_exon_length > longest_transcript: 
        longest_transcript = current_exon_length
        longest_transcript_ID = old_transcript_ID
    current_exon_length = 0
    return(current_exon_length,longest_transcript,longest_transcript_ID)


        
def exon_length_filter(file_name = "test",source_pathway_name = os.getcwd(),deposit_pathway_name =os.getcwd(),gen_dict = {"ENSG00000160072":["ENST00000673477","ENST00000472194","ENST00000378736","ENST00000308647","ENST00000442483"],"ENSG00000225972":["ENST00000416931"],"ENSG00000279928":["ENST00000624431","ENST00000424215"],"ENSG00000142611":["ENST00000378391","ENST00000607632","ENST00000511072"]}):
    """This function selects only the transcripts that have the longest total mRNA"""  
    
    print("Representative transcipts are filtered based on exon length. Please wait...")
    bar,start_time = te.bar_builder(length_multiplyer = 3)
    source_pathway_name,deposit_pathway_name = te.__do_pathways_exist__(source_pathway_name,deposit_pathway_name)
    total_genes = len(gen_dict)
    gens_done = 0

    with open(source_pathway_name+"\\"+file_name+".gtf", 'r') as f:
        
        old_gen = str()
        old_transcript_ID = str()
        representative_transcript = dict()
        representative_trasnscript_not_found = True
        longest_transcript_ID = str()
        current_exon_length = 0
        longest_transcript = 0 
        percentage_done = 0
        
        for entry in f: 
            
            try:
                corrent_gen = te.gene_ID_finder(entry)
            except:
                corrent_gen = old_gen
            #The block above test if there is a gen name in the entry
            if corrent_gen != old_gen:   
                representative_trasnscript_not_found = True

            #The block above determines if the Gen name is new and set the test
            #representative_trasnscript_not_found back to true which is used to 
            #make the program faster if there is just one transcript for a given
            #gen in the dict
            if representative_trasnscript_not_found and corrent_gen != str():
                #print(corrent_gen)
                #The conditon prvents serges if a representative transcript has
                #all ready been chosen
                if corrent_gen != old_gen:
                    current_exon_length,longest_transcript,longest_transcript_ID = __longest_transcript_finder(current_exon_length,longest_transcript,longest_transcript_ID,old_transcript_ID)
                    representative_transcript[old_gen] = longest_transcript_ID
                    try:
                        del gen_dict[old_gen]
                        old_gen = corrent_gen                   
                        gens_done += 1
                        corrent_percentage_done = (gens_done/total_genes)*100
                        if corrent_percentage_done > percentage_done+10:
                            bar,start_time = te.bar_builder(percentage=percentage_done+10,length_multiplyer = 3,start_time=start_time,bar =bar)
                            percentage_done = int(corrent_percentage_done)  
                        
                         
                    except:
                        old_gen = corrent_gen
                    longest_transcript = 0
                    #The block above adds the transcript of the last gen that 
                    #had the longest exons into the representative transcripts dict
                    try: 
                        #This try / except block test if the gen is in the input dictionary
                        transcript_IDs = gen_dict[corrent_gen]
                        if len(gen_dict[corrent_gen]) == 1:
                            #This conditions is a short cut for Genes that 
                            #allready have a representative transcript
                            representative_transcript=gen_dict[corrent_gen[0]]
                            representative_trasnscript_not_found = False
                            continue
                    except:
                        continue
                    
                try: 
                    current_transcript_ID = te.transcript_ID_finder(entry)         
                except: 
                    continue
                #The block above searches for a transcript ID in the  current enty

                if current_transcript_ID in transcript_IDs:
                    #This condition test if the Transcript is one of the 
                    #candidates for representative transcripts
                    if current_transcript_ID != old_transcript_ID:
                        #This condition if the enty still belongs to the 
                        #previous transcript and is triggers if that is not the case
                        current_exon_length,longest_transcript,longest_transcript_ID = __longest_transcript_finder(current_exon_length,longest_transcript,longest_transcript_ID,old_transcript_ID)
                        try:
                            transcript_IDs.remove(old_transcript_ID)
                            old_transcript_ID = current_transcript_ID
                        except:
                            old_transcript_ID = current_transcript_ID
                    if exon_fider(entry): 
                        exon_length = exon_length_calculator(entry)
                        current_exon_length += exon_length
                    else: 
                        continue 
        current_exon_length,longest_transcript,longest_transcript_ID = __longest_transcript_finder(current_exon_length,longest_transcript,longest_transcript_ID,old_transcript_ID)
        representative_transcript[old_gen] = longest_transcript_ID
    del representative_transcript[str()]
    te.bar_builder(100,length_multiplyer = 3,start_time=start_time,bar =bar)
    return(representative_transcript)

if __name__ == "__main__":
    exon_length_filter()
    
