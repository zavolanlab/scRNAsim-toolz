# cDNA Generator module

Description of the module:
The function of this module is to generate cdDNA based on mRNA transcript seuqences and the coresponding priming probabilities. 


**Input files**


transcript_copies (csv-formatted) containing:

- ID of transcript
- ID of parent transcript
- copy number


_Eample_

`[ID of transcript]    [ID of parent transcript]    [copy number]`


transcript_sequences (fasta-formatted) containing:
 
- ID of transcript 

_Eample_

`> [ID of transcript]
AGUGACGUUAGACCAGAUAGAC....`


priming_sites (gtf-formatted) containing:

- ID of transcript
- Position of priming site
- Binding likelihood of priming site

_Eample_

`[ID of transcript]    ... [Position of priming site]... [Binding likelihood of priming site]`


**Output files**

cDNA (fasta-formatted) containing:

- cDNA sequence ID
- Uniquie cDNA sequence



cDNA_counts (csv-formatted) containing:

- cDNA sequence ID
- cDNA copy number






