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

import pandas as pd
import math

class CreateTranscript():
    def __init__(self):
        with open("RIBlast output example.txt", 'r') as file:
            self.raw_interactions = file.readlines()[3:]

    
    
    def generate_interaction_df(self):
        
        self.interaction_list = []
        
        #clean up the original list so that we have a neet list
        for i in range(0, (len(self.raw_interactions)-1)):
            current_interaction = self.raw_interactions[i].strip(' \n').replace('(', '').replace(')','').replace('-',',').replace(':',',').split(',')
            self.interaction_list.append(current_interaction) 
        
        #identify if the interaction is the same as the previous one, just shifted by 1 bp (if we have 20 A in the transcript the 15 T primer has 5 matching possibilities although it is only 1 bindingsite)
        previous_interaction_base = int
        for i in range(0, len(self.interaction_list)):
            previous_interaction_base = int(self.interaction_list[i-1][13])
            if int(self.interaction_list[i][13]) in range(previous_interaction_base-1,previous_interaction_base-15,-1):
                self.interaction_list[i].append('Repeat')
            else :
                self.interaction_list[i].append('Not_repeat')

        #exclude all interactions which are a repeat and belong to the same bindingsite
        self.cleaned_interaction_list = [item for item in self.interaction_list if item[-1]=='Not_repeat']
        
        
        #add total number of interactions per transcript and calculate energy
        self.df = pd.DataFrame(self.cleaned_interaction_list)
        self.df['Number of interactions'] = int
        self.df['Interaction Energy'] = float
        
        energy_constant = 1.380649*10**(-23)*298
        kcalmol_joul = 6.9477*10**-21
        
        
        for ind in self.df.index:
            self.df['Number of interactions'][ind]=self.df[3].value_counts()[self.df[3][ind]]
            self.df['Interaction Energy'][ind]=math.exp(-float(self.df[5][ind])*kcalmol_joul/energy_constant)
        print(self.df['Interaction Energy'])
        print(self.df)
        

        return self.df



transcripts = CreateTranscript()    
interaction_df = transcripts.generate_interaction_df()        

#print line by line to file and then you're done

output = str()
for i in interaction_df.index:
    #print(interaction_df[3][i]+'\t' + 'RIBlast' + '\t' + 'Priming_site' + '\t' + interaction_df[13][i] + '\t' + interaction_df[12][i] + '\t' + '.' + '\t' + '+' + '\t' + '.' + '\t' + f'Accessibility_Energy "{interaction_df["Interaction Energy"][i]}"')
    output = output + str(interaction_df[3][i]+'\t' + 'RIBlast' + '\t' + 'Priming_site' + '\t' + interaction_df[13][i] + '\t' + interaction_df[12][i] + '\t' + '.' + '\t' + '+' + '\t' + '.' + '\t' + f'Accessibility_Energy "{interaction_df["Interaction Energy"][i]}"' + '\n')

print(output)
with open('output_transcripts_df.txt', 'w') as f:
    f.write(output)

