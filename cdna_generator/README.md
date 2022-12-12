# cDNA Generator module

Description of the module:
The function of this module is to generate cdDNA based on mRNA transcript seuqences and the coresponding priming probabilities. 

**Example usage**

    python ../cdna/cli.py -ifa yeast_example.fa -icpn copy_number_input.csv -igt Example_GTF_Input.GTF -ofa cDNA.fasta -ocsv cDNA.csv


**Input files**


transcript_copies (csv-formatted) containing:

- ID of transcript
- ID of parent transcript
- transcript copy number


transcript_sequences (fasta-formatted) containing:
 
- ID of transcript 
- transcript-sequence

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






