# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:48:04 2022

@author: baerma
"""

# with open("RIBlast output example.txt", 'r') as file:
#     content = file.readlines()
#     print(content[3:][0].strip(' \n').split(',')[-1].strip('()').split(':'))
#     #create a instant of each transcript class
#     print(len(content[3:]))
#     print((content[3:]))

class CreateTranscript():
    def __init__(self):
        with open("RIBlast output example.txt", 'r') as file:
            self.list_of_interactions = file.readlines()[3:]

    
    
    def generate_interaction_list(self):
        interaction_list = []
        for i in range(0, (len(self.list_of_interactions)-1)):
            current_interaction = self.list_of_interactions[i].strip(' \n').split(',')
            #print(self.list_of_interactions[i].strip(' \n').split(','))
            interaction_list.append(current_interaction)
        return interaction_list

            

transcriptlist = CreateTranscript()            
print(transcriptlist.generate_interaction_list())

#go from interaction list to transcript list? -no we will serve them a interaction list.
            
  #          print(content[3:][0].strip(' \n').split(',')[-1].strip('()').split(':'))

