'''
This function writes keeps the representative transcripts from the original input file (gtf)
and writes them to an output.
(The representative transcripts being listed in a csv file)
'''

def gtf_representative_transcripts(original_gtf, transcripts):
    transcript = pd.read_csv(transcripts)
    representative_transcripts = []

    with open (original_gtf, 'r') as file:
        for id in transcript['id']:
            for line in file:
                if id in line:
                    representative_transcripts.append(line)


    with open ('output.gtf', 'w') as outputfile:
        outputfile.write(representative_transcripts)


