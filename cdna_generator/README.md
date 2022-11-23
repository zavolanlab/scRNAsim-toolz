# cDNA Generator module

Description of the module:
The function of this module is to generate cdDNA based on mRNA transcript seuqences and the coresponding priming probabilities. 


**Input files**


transcript_copies (csv-formatted) containing:

- ID of transcript
- ID of parent transcript
- copy number


_Eample_

`[id of generated transcript]    [ID]    [Count]`


transcript_sequences (fasta-formatted) containing:
 
- id of generated transcript? (in the header)

_Eample_

`> [id of generated transcript]
AGUGACGUUAGACCAGAUAGAC....`


priming_sites (gtf-formatted) containing:

- id of generated transcript?
- position of priming site and binding likelihood 

_Eample_

`[id of generated transcript]    ... [position of priming site]... [binding likelihood ]`


**Output files**

cDNA (fasta-formatted) containing:

- cDNA sequence ID
- uniquie cDNA sequence and "cDNA sequence ID



cDNA_counts (csv-formatted) containing:

- cDNA sequence ID
- cDNA count






