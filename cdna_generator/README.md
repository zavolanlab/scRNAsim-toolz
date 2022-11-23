# cDNA Generator module

Description of the module:
The function of this module is to generate cdDNA based on mRNA transcript seuqences and the coresponding priming probabilities. 


**Input files**


transcript_copies (csv-formatted) containing:

- ID of transcript
- ID of parent transcript
- copy number


transcript_sequences (fasta-formatted) containing:
 
- ID of transcript 


priming_sites (gtf-formatted) containing:

- ID of transcript
- Position of priming site
- Binding likelihood of priming site



**Output files**

cDNA_sequences (fasta-formatted) containing:

- cDNA sequence ID
- cDNA-sequence


cDNA_counts (csv-formatted) containing:

- cDNA sequence ID
- cDNA-counts






