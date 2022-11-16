# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 16:47:20 2022

@author: baerma
"""


with open("RIBlast output example.txt", 'r') as file:
    content = file.readlines()
    print(content[3:][0].strip(' \n').split(',')[-1].strip('()').split(':'))
    #create a instant of each transcript class




#os.chdir('C:/Users/baerma/Desktop/PhD-Local/Lectures/Programming for Life Sciences/priming-site-predictor/primingsitepredictor')

