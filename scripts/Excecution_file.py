### Imports ###
import os

import transkript_extractor as te
import Exon_length_filter as elf
import representative_v4 as rtcl

### Scipt ###
def exe(file_name = "test",source_pathway_name = os.getcwd(),deposit_pathway_name = os.getcwd(),Input_free = True):
    file_name,source_pathway_name_2,deposit_pathway_name_2 = te.extract_transkript(file_name,source_pathway_name,deposit_pathway_name,Input_free = Input_free)
    inter_mediate_file_directory = os.path.join(deposit_pathway_name,file_name+"_intermediate_file.txt")
    print("Transcripts are filterd based on transcipt score please wait...")
    pre_filter_representative_transcripts_dict = rtcl.find_repr_by_SupportLevel(inter_mediate_file_directory)
    print("Transcripts filtered\n")
    elf.exon_length_filter(file_name,source_pathway_name,deposit_pathway_name,gen_dict= pre_filter_representative_transcripts_dict,Input_free = Input_free)
    return(file_name,source_pathway_name,deposit_pathway_name)
### from consol ####
##D:\\Uni\\Sem 9\\Programing in the Life sciences\\Projekt\\Intermediat Files
if __name__ == "__main__":
    exe()