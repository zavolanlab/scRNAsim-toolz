#### Transcript extractor #####

### Called Packages ###
import re
import os

### Functions ###
    
def __searche_for_preexisting_files(file_name,deposit_pathway_name = os.getcwd()):
    """This function searches for preexisting files of the same name as the results file of the current program. It allows the user to choose to move on with the pre-existing file """
    File_of_same_name_found = False
    generat_new_file = False
    directory_content = os.listdir(deposit_pathway_name)
    for file in directory_content: 
        Search_profile = file_name+"_intermediat_file.txt"
        if file == Search_profile: 
            while True: 
                File_found_input = input ("An intermediate file has allready been generated from this file\nDo you want to generate a new one [y/n] \n>")
                if File_found_input == "n":                     
                    File_of_same_name_found = True
                    break
                elif File_found_input == "y":
                    generat_new_file = True
                    break
                else: 
                    print("Sorry this was not a valid input\nPlease press [y] if you want to generate a new file or [n] if you want to use the preexisting file")
            break
        else: 
            continue
    if File_of_same_name_found: 
        print("No new file will be generated, the program can continue")
    elif generat_new_file: 
        print("A new file will be generated please wait...\n")
    else:            
        print("No pre-existing intermediate file based on the currend file have been found.\nA new file will be generated please wait...\n")
    return(File_of_same_name_found)

def __test_file_name(file_name,source_pathway_name = os.getcwd()):
    """This function validates that the source file exists at the source path. It turns the file name input in a standardized format that can be used in the next steps"""
    
    directory_content = os.listdir(source_pathway_name)
    
    index_of_the_dot = file_name.rfind(".")
    valide_source_file = False
    validate_source_file = True
    if index_of_the_dot ==-1:
        file_name += ".gtf"       
    else: 
        source_file_typ = file_name[index_of_the_dot:]
        not_a_file_type = re.compile(".\d{1,13}")
        try_not_a_file_type = not_a_file_type.search(source_file_typ)
        if source_file_typ == ".gtf":
            file_name = file_name
        elif try_not_a_file_type:
            file_name += ".gtf"
        else: 
            print("This program can not handle",source_file_typ,"files. \nplease use a .gtf file" )
            validate_source_file = False
    #The block above tests if the file_name includes the file type and if no 
    #file type is found adds ".gtf" und if a non ".gtf" file is found gives an error
    
    if validate_source_file: 
        for file in directory_content: 
            if file == file_name:
                valide_source_file = True 
                break
    #The block above tests if a file on the given name is in the given directora 
    
    if valide_source_file:
        print("The file:",file_name,"has been found.\n")
    else: 
        print("No .gtf file of the name",file_name,"has been found in this pathway")
    #The bock above gives feed back regarding the results of the file test 
    
    file_name = file_name.replace(".gtf","")
    #This line normalizes the file name 
    return(valide_source_file,file_name)

def __do_pathways_exist__(source_pathway_name,deposit_pathway_name):
    """This funtion tests that the entered pathways actualy exist"""
    does_source_pathway_exist = os.path.exists(source_pathway_name)
    does_deposit_pathway_exist = os.path.exists(deposit_pathway_name)
    #The Block above does the actual testing
    if does_source_pathway_exist:
        source_pathway_name = source_pathway_name
    else: 
        print("The source pathway:",source_pathway_name,"has not been found\nThe source pathway was set to the default")
        source_pathway_name = os.getcwd()
    #The block above detail the possible reactions for the source pathe existing or not existing
    if does_deposit_pathway_exist: 
        deposit_pathway_name = deposit_pathway_name
    else: 
        print("The deposit pathway:",deposit_pathway_name,"has not been found\nThe deposit pathway was set to the default")
        deposit_pathway_name = source_pathway_name
    #The block above details the possible reactions for the deposit pathway existing or not existing 
    return(source_pathway_name,deposit_pathway_name)
        
def gene_ID_finder(entry):
    """This function is supposed to finde the gen ID of a known gen entry"""
    index_gene_id = entry.find("gene_id")
    #This line determines where the transcript ID is 
    find_gene_id_name = re.compile("\"\S{1,25}\"")
    #This line defines the pattern of a transcript ID
    sub_entry = entry[index_gene_id:]
    #This line generates a subentra starting at the locatind of the tanscript ID
    try_find_gene_id_name = find_gene_id_name.search(sub_entry)   
    gene_ID = try_find_gene_id_name[0].replace("\"","")
    #The block above findes the transcript ID and changes it into a usable format
    return (gene_ID)
       
def transcript_ID_finder (entry):
    """This function is supposed to finde the transcript ID in a known transcript entry"""
    index_transcript_id = entry.find("transcript_id")
    #This line determines where the transcript ID is 
    find_transcript_id_name = re.compile("\"\S{1,25}\"")
    #This line defines the pattern of a transcript ID
    sub_entry = entry[index_transcript_id:]
    #This line generates a subentra starting at the locatind of the tanscript ID
    try_find_transcript_id_name = find_transcript_id_name.search(sub_entry)   
    
    try: 
        transcript_ID = try_find_transcript_id_name[0].replace("\"","")
        #The block above findes the transcript ID and changes it into a usable format
    except:
        transcript_ID = ""
    return (transcript_ID)
        
def transcript_support_level_finder(entry):
    """This function is supposed to find the transcript support level in a known transcript entry"""
    transcript_support_level_start_ID = entry.find("transcript_support_level")
    sub_entry = entry[transcript_support_level_start_ID:]
    
    try:
        #the try and except clauses are there for the case that there is 
        #no actual score given
        score_fidner = re.compile("\W\w{1,16}\W{2}")
        try_score_finder = score_fidner.search(sub_entry)              
        Pre_score_1 = try_score_finder[0]
        Pre_score_2 = Pre_score_1.replace("\"","")
        Pre_score_2 = Pre_score_2.replace("(","")
        transcript_support_level = Pre_score_2.replace(";","")
        if transcript_support_level == "NA":
            transcript_support_level = 100

    except:
        #will give NA if no numerical score is found 
        transcript_support_level = 100
    return (transcript_support_level)



    
def _transcript_extractor (file_name,source_pathway_name,deposit_pathway_name): 
    """This funtion extracts the transcript number ,transcript ID, the transcript support level, the transcrip length and the line index from a gtf file of a given name and saves tham as a new file name given_name_intermediat_file.txt. It only works in the directory of the skript at this point"""
    with open(source_pathway_name+"\\"+file_name+".gtf", 'r') as f:
        #this line opens the file in question 
        
        Old_gen_ID = str() 
        #stand-in as the first couple entrys are not genes
        with open(deposit_pathway_name+"\\"+file_name+"_"+"intermediat_file"+".txt","w") as IMF:
            #creates a new file that will have the same name as the source file 
            #but with an added _inter_mediat_file
            transcript_number = 0
            for entry in f: 
                #this loop reads all lines in the source file one by one 
                Gen_finder = re.compile("gene_name") 
                try_gen_finder = Gen_finder.search(entry)
                #the lines above determin if the is a "gene_name" collumn 
                #in the current entry
                if (try_gen_finder):
                    Gen_ID = gene_ID_finder(entry)
                else:
                    Gen_ID = Old_gen_ID
                #determins the gene_name using the function showen futher up 
                if Gen_ID != Old_gen_ID:
                    #chachs if the gene is new or if it is an entry belonging 
                    #to the previous gen e.g. a transkript, exon ...
                    Gen_entry = ">"+str(Gen_ID)+"\n"
                    IMF.write(Gen_entry)
                    #if the Gen is new a new line is inserted in the intermediat 
                    #file to seperate the transcripts belonging to the different genes
                    transcript_number = 0
                    #the transcript number indicates how many transcripts have
                    #been found of the same gene befor the current one it starts 
                    #at 1 for each new gene 
                    Old_gen_ID = Gen_ID
                    #the new_name is set as the Old meaning current name 
                
                transkript_finder = re.compile("\ttranscript\t")
                try_transkript_finder = transkript_finder.search(entry)
                #the lines above serve to identify if the current entry blonges 
                #to a transcript 
                
                if try_transkript_finder: 
                    #this confition is activated if the entry belongs to a transcript 
                    transcript_number += 1
                    #the transcript number is updated (starting at on given that t
                    #he default is 0)
                    Transcript_ID  = transcript_ID_finder (entry)
                    #the function that determins the transcript ID is called
                    transcript_support_level = transcript_support_level_finder(entry)
                    #the function that determins the transcript support level is called
                    New_entry = str(transcript_number)+"\t"+str(Transcript_ID)+"\t"+str(transcript_support_level)+"\t"+"\t\n"
                    IMF.write(New_entry)
                    #this lines assemble the transcript number ,transcript ID 
                    #and the transcript support level into a new line for the 
                    #intermediat file and adds them in
        print("The transcripts have been collected")
        
def __parameter_editor(file_name,source_pathway_name,deposit_pathway_name):
    """This funtion allows for chaging the parameters after running the program"""
    while True:
        #This is the main loot wich allows show the current partameters and 
        #allows the usert to shouse if the want to use edit the paremeters or 
        #to continue with the given parameters
        print("The program will run with the following parameters:\nFile name:\t\t",file_name,"\nSource pathway:\t",source_pathway_name,"\nDeposit pathway:\t",deposit_pathway_name,"\n")
        parameter_conformation = input("To continue with these parameters input [continue or c] to change them input [edit]\n>")
        #The block above incluts  the "table" whith the current parameters as well as the input prompt
        if parameter_conformation == "continue"or parameter_conformation =="c":
            #This if clause will end the main loop and the funtion 
            break
        elif parameter_conformation == "edit":
            #This conditin lead into the secundary loop which allows the user to edit the parameters
            while True: 
                change_question = input("select the parameter you want to change [nfile/spath/dpath] or input [b] to go back\n>")
                #This is the input prompt for the secundary loop 
                if change_question == "nfile":
                    #This condition allows the user to chenge the file name 
                    file_name = input("Pleas input the new file name\n>")
                    break
                elif  change_question == "spath":
                    #This condition allows the user to change the source path
                    source_pathway_name = input("Pleas input the new source path\n>")
                    
                    does_source_pathway_exist = os.path.exists(source_pathway_name)
                    if does_source_pathway_exist:
                        break
                    else: 
                        print("The new source pathway:",source_pathway_name,"does not existe\nThe source pathway was returnt to default:",os.getcwd())
                        source_pathway_name = os.getcwd()
                    #The block above test that if the new source pathway is valid
                elif  change_question == "dpath":
                    #This condition allows the user to change the deposit path
                    deposit_pathway_name = input("Pleas input the new deposit path name\n>")
                    
                    does_deposit_pathway_exist = os.path.exists(deposit_pathway_name)
                    if does_deposit_pathway_exist:
                        break
                    else:
                        print("The new deposit pathway:",deposit_pathway_name,"does not existe\nThe deposit pathway was returnt to default:",source_pathway_name)
                        deposit_pathway_name = source_pathway_name
                    #The block above test if the new deposit pathway is valid
                elif  change_question == "b":
                    # This condition allows the user to return to the main loop
                    break             
                else:
                    #This condition covers all non valid inputs into the secund loop
                    print("The input",change_question,"is not valide please use one of the specifyed comands") 
                    
        else: 
            #This condition covers all non valid input for the main loop 
           print("The input",parameter_conformation,"is not valide please use one of the specifyed comands\n") 
    return(file_name,source_pathway_name,deposit_pathway_name)
        
def extract_transkript (file_name = "test",source_pathway_name = os.getcwd(),deposit_pathway_name = False): 
   """This it the overall exetutable funtion that will execute the transcript extraction process for a given file with all checks. 
   The default file name is "test". This function will also return the file name, the source pathway and the depisti pathway that have been used to generate the intermediat file"""
   if deposit_pathway_name == False: 
       deposit_pathway_name = source_pathway_name  
   file_name,source_pathway_name,deposit_pathway_name = __parameter_editor(file_name,source_pathway_name,deposit_pathway_name)
   source_pathway_name,deposit_pathway_name =__do_pathways_exist__(source_pathway_name,deposit_pathway_name)
   validated_file_name = __test_file_name(file_name,source_pathway_name)
   file_name = validated_file_name[1]
   if validated_file_name[0]:
       if __searche_for_preexisting_files(file_name,deposit_pathway_name):
           print("The transcripts has been collected\n")
       else:
           _transcript_extractor (file_name,source_pathway_name,deposit_pathway_name)
   return(file_name,source_pathway_name,deposit_pathway_name)

#### Dev part ####

if __name__ == "__main__":
    extract_transkript()
#This line allows the file to be executed on its own also from 


