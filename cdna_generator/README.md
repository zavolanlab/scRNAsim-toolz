# cDNA Generator module

Description of the module:
The function of this module is to generate cdDNA based on mRNA transcript seuqences and the coresponding priming probabilities. 


**Input files**

Copy_number_file:

- csv-formatted file ("NewTranscriptID,ID,Count")

- id of generated transcript

- id of original transcript (without intron inclusions)
count

_Eample_

`[id of generated transcript]    [ID]    [Count]`


transcript_sequences_file:

- fasta-formatted file 

- id of generated transcript? (in the header)

_Eample_

`> [id of generated transcript]
AGUGACGUUAGACCAGAUAGAC....`


priming_site_file:

- gtf-formatted file 

- id of generated transcript?

- position of priming site and binding likelihood 

_Eample_

`[id of generated transcript]    ... [position of priming site]... [binding likelihood ]`


**Output files**

cDNA_file:

- fasta-formatted file 

- Includes all the uniquie "cDNA sequence" and "cDNA sequence ID"



cDNA_count_file:

- csv-formatted file 

- Includes "cDNA sequence ID" and "cDNA count"






