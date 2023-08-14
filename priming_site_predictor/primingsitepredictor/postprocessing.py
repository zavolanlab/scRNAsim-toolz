# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 14:06:20 2022

@author: baerma
"""

import pandas as pd
import math
pd.options.mode.chained_assignment = None

class PostProcessRIBlast():
    
    def __init__(self):
        output = self.generate_gtf()
        #print(output)
     
    
    def calculate_energy(self, value):
        energy_constant = 1.380649*10**(-23)*298
        kcalmol_joul = 6.9477*10**-21
        return (math.exp(-float(value)*kcalmol_joul/energy_constant))
    
    
    def create_list_from_output(self):
        self.file = "RIBlast output example.txt"
        self.firstline = 3
        self.interaction_list = []
        with open(self.file, 'r') as file:
            self.raw_interactions = file.readlines()[self.firstline:]   
        self.number_entries = len(self.raw_interactions)
        
        for i in range(0, self.number_entries-1):
            current_interaction = self.raw_interactions[i].strip(' \n').replace('(', '').replace(')','').replace('-',',').replace(':',',').split(',')
            self.interaction_list.append(current_interaction) 

        return self.interaction_list
            

    def create_pandas_df(self):
        self.interaction_list = self.create_list_from_output()
        self.df = pd.DataFrame(self.interaction_list)
        self.df['Number_of_interactions'] = int(0)
        self.df['Interaction_Energy'] = float(0)
        self.transcript = 3
        self.energy = 5
        
        for index in self.df.index:
            self.df['Number_of_interactions'][index]=self.df[self.transcript].value_counts()[self.df[self.transcript][index]]
            self.df['Interaction_Energy'][index]=self.calculate_energy(self.df[self.energy][index])

        self.df['Normalised_interaction_energy']=self.df['Interaction_Energy']/self.df['Number_of_interactions']
        
        return self.df

    
    def generate_gtf(self):
        self.interaction_df = self.create_pandas_df()
        self.output = str()

        for index in self.interaction_df.index:
            self.output = self.output + str(self.interaction_df[3][index]+'\t' + 'RIBlast' + '\t' + 'Priming_site' + '\t' + self.interaction_df[13][index] + '\t' + self.interaction_df[12][index] + '\t' + '.' + '\t' + '+' + '\t' + '.' + '\t' + f'Interaction_Energy "{self.interaction_df["Normalised_interaction_energy"][index]}"' + '\n')
        
        with open('output_transcripts_df.GTF', 'w') as f:
            f.write(self.output)
            return(self.output)
        

print(PostProcessRIBlast().output)