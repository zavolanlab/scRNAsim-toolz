#### Find representative transcripts ####
"""Version 1.1.1"""

### Imports ### 
import argparse

### Functions ###
    
def attributs_converter(attributs):
    """
    This funtion converts the "unstrucktured" ;-seperated part of he line into a list of identifyers and coresponding data the struckture of
    which can be used ot find the data easyly e.g the index of the identifier transcrip_id + 1 will give the trasncript id of the current gene
    Input: 
        attributs = str() #the unstrucktured part of the entry
    Output:
        attributs = list() # cleand list with the characterritsics discribed above
    """
    attributs = attributs.replace("\"","")
    attributs = attributs.replace(";","")
    attributs = attributs.replace("\\n","")
    attributs =attributs.split(" ")
    
    return(attributs)

def find_in_attributs (attributs,look_for):
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
        index = attributs.index(look_for)+1
        return attributs[index]
    except: 
        #print("No",look_for,"in the entry the return was set to NA\n",attributs)
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
        
    

def get_rep_trans(file_name = "test"):
    """ 
    This is the main function of this script it selects one representative transcrip per gene based on a gtf annotation file. 
    It does so be two criteria: first the transcript support level and it there are several transcript 
    of one gene that have the same trasncript_support_level it chooses the one that corresponds to the longest mRNA.
    Input: 
        file_name = str() # name of the annotation file with or without the .gtf part
    Output: 
        rep_transcripts = {gene_id : transcript_id}
    """
    
    #setting defoult variables
    rep_trans = dict()
    cur_gID = str()
    cur_best_trans = [str(),100,0] # [transcript_id , transcript_support_level , transcript_length]
    pot_best_trans = False
    cur_tID = str()
    ignor_trans = False
    
    with open (file_name,"r") as f: 
        for line in f: 
            entry = line.split("\t")
            
            #removes expected but unneeded entrys
            exp_unneed = ["CDS","stop_codon","five_prime_utr","three_prime_utr","start_codon",'Selenocysteine']
            if len(entry) == 1 or entry[2] in exp_unneed:
                continue
            
            #this function turns the less organized part of the entry into a reable list
            attributs = attributs_converter(entry[8])
            #looking for and processing exons entrys
            if entry[2] == "exon": 
                
                #dicide if to contiune or not
                if ignor_trans: 
                    continue
                elif cur_gID != attributs[1]:
                    raise ValueError("ERROR exon from an unexpected Gen")
                    continue
                elif find_in_attributs (attributs,"transcript_id") != cur_tID:
                    raise ValueError("exon from an unexpected transcript")
                    continue
                
                #adding the length of the exon to the appropriat list and chacking for changes in best transcript
                if pot_best_trans: 
                    pot_best_trans[2]+= int(entry[4])-int(entry[3])
                    if pot_best_trans[2] > cur_best_trans[2]: 
                        cur_best_trans = pot_best_trans
                        pot_best_trans = False
                else:
                    cur_best_trans[2]+= int(entry[4])-int(entry[3])

                                       
                
            #looking for and processing transcript entrys
            elif entry[2] == "transcript":
                    
                #varryfi that the gen is correct
                if cur_gID != attributs[1]:
                    raise ValueError("ERROR transcript from an unexpected Gen")
                    continue
                
                #finding the transcript id and the support level
                cur_tID = find_in_attributs (attributs,"transcript_id")       
                t_supp_lvl = find_in_attributs (attributs,"transcript_support_level")    
                
                #If there is no transcript support level or the level is given as NA it is nomed as 100. else the transcript support level is tunrn into int
                if t_supp_lvl == "NA": 
                    t_supp_lvl = 100
                else:
                    try:
                        t_supp_lvl = int(t_supp_lvl)
                    except: 
                        t_supp_lvl = 100
                
                
                #decides if the transcript has potential to become the representative transcript
                if t_supp_lvl < cur_best_trans[1] or cur_best_trans[0] == "":
                    cur_best_trans = [cur_tID,t_supp_lvl,0]
                    pot_best_trans = False
                    ignor_trans = False
                     
                elif t_supp_lvl == cur_best_trans[1]:
                    pot_best_trans = [cur_tID,t_supp_lvl,0] 
                else:
                    ignor_trans = True
                
                  
            #looking for and processing gene entrys
            elif entry[2] == "gene":
                
                #updating rep_trans dict
                if cur_gID not in rep_trans: 
                    rep_trans[cur_gID] = cur_best_trans
                else: 
                    if rep_trans[cur_gID][1] > cur_best_trans[1]: 
                        rep_trans[cur_gID] = cur_best_trans
                    elif rep_trans[cur_gID][1] == cur_best_trans[1] and rep_trans[cur_gID][2] < cur_best_trans[2]: 
                        rep_trans[cur_gID] = cur_best_trans
                
                #updating cur_gID and resetting cur_best_trans
                cur_gID = attributs[1]
                cur_best_trans = [str(),100,0]
                    
            #raises an error for unidentifyable entrys
            else: 
                raise ValueError("This entry could not be identified\n",entry)
        
        #addding the final gene to the dictionary
        if cur_gID not in rep_trans: 
            rep_trans[cur_gID] = cur_best_trans
        else: 
            if rep_trans[cur_gID][1] > cur_best_trans[1]: 
                rep_trans[cur_gID] = cur_best_trans
            elif rep_trans[cur_gID][1] == cur_best_trans[1] and rep_trans[cur_gID][2] < cur_best_trans[2]: 
                rep_trans[cur_gID] = cur_best_trans        
        
        del rep_trans[""]
        rep_transcripts = _re_format(rep_trans)
        return(rep_transcripts )

def gtf_file_writer (original_file, output_file): 
    """
    this function writes the output GTF file
    """
    output = []
    rep_transcript_dict = get_rep_trans(original_file)

    with open(original_file, 'r') as f:
            for entry in f: 
                if entry[0] != '#':
                    attributes = attributs_converter(entry)
                    type_ = attributes[2]
                    if type_ == 'gene':
                        gene_id = find_in_attributs(attributes, 'gene_id')
                        output.append(entry)
                    if type_ != 'gene':
                        transcript_id = find_in_attributs(attributes, 'transcript_id')
                        if rep_transcript_dict[gene_id] == transcript_id:
                            output.append(entry)

    with open(output_file, 'w') as last_file:
        last_file.write(output)

def _test(): 
    """
    This funtion is ment to be run for test
    Output: 
        file with the dictionary generated based on the test file 
    """
    file_name = "test.gtf"
    rt = get_rep_trans(file_name)
    expected_result = {"ENSG00000160072":"ENST00000472194","ENSG00000234396":"ENST00000442483",
                       "ENSG00000225972":"ENST00000416931","ENSG00000224315":"ENST00000428803",
                       "ENSG00000198744":"ENST00000416718","ENSG00000279928":"ENST00000624431",
                       "ENSG00000228037":"ENST00000424215",'ENSG00000142611':'ENST00000378391'}
    if rt != expected_result: 
        print("The test fail due to not yieding the same results")
        print("The results the program got\n",rt)
        print("The expected results\n",expected_result)
    else: 
        print("The test was succses full")
            
### Execution part ###
if __name__ == "__main__":   
    parser = argparse.ArgumentParser(description="find_representativ_transcripts",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-file_name", required=True, help="gtf file with genome annotation")
    parser.add_argument("-t", required=False,default = False,help="to run the test input -t True")
    args = parser.parse_args()
    
    #standadize the file_name inlude .gtf#

    file_name = args.file_name
    i_gtf = file_name.find(".gtf")
    if i_gtf == -1:
        file_name += ".gtf"  
    
    if args.t: 
        _test()
    else:
        get_rep_trans(file_name)
   

    
                
            