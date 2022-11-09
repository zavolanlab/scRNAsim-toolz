import pandas as pd
import numpy as np

'''
Sample transcript 

This part of the code does Poisson sampling proportionally to gene expression levels for each gene. 
 
input:  total transcript number (int) 
        csv file with gene id and  gene expression levels (columns named 'id' and 'level')

output: csv file with gene id and count
        gtf file with transcript samples
'''

def transcript_sampling(total_transcript, transcripts):

    #read file containing representative transcript levels and id

    transcript = pd.read_csv(transcripts)
    levels = []

    #poisson sampling for each gene, proportional to expression levels

    for expression_level in transcript['level']:

        poisson_sampled = np.random.poisson(total_transcript/expression_level)
        levels.append(poisson_sampled)
        # note: if levels from input are decimals, total trancript should be multiplied by expr level, not divided
        # poisson_sampled = np.random.poisson(total_transcript*expression_level)

    #write output csv file containing transcript id and count (representative transcript numbers)

    transcript_numbers = pd.DataFrame({'id': transcript['id'],'count': levels})
    pd.DataFrame.to_csv(transcript_numbers, "representative_transcript_numbers.csv")


