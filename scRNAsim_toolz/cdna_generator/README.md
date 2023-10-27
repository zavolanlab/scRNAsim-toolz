# cDNA Generator module
## Usage
```
usage: cdna-generator [-h] -ifa INPUT_FASTA -igtf INPUT_GTF -icpn INPUT_COPY_NUMBER -ofa OUTPUT_FASTA -ocsv OUTPUT_CSV [-v]

Generate cDNA sequences based on primer probabilities.

options:
  -h, --help            show this help message and exit
  -ifa INPUT_FASTA, --input_fasta INPUT_FASTA
                        genome fasta file
  -igtf INPUT_GTF, --input_gtf INPUT_GTF
                        gtf file
  -icpn INPUT_COPY_NUMBER, --input_copy_number INPUT_COPY_NUMBER
                        input copy number (csv) file
  -ofa OUTPUT_FASTA, --output_fasta OUTPUT_FASTA
                        output fasta file
  -ocsv OUTPUT_CSV, --output_csv OUTPUT_CSV
                        output fasta file
```
Example:
```
cdna-generator -ifa tests/cdna_generator/files/transcript.fasta -igtf tests/cdna_generator/files/Example_GTF_Input.GTF -icpn tests/cdna_generator/files/copy_number_input.csv -ofa cdna_seq.fa -ocsv cdna_counts.csv
```

## Overview
Generate cDNA based on mRNA transcript sequences and the coresponding priming probabilities. 


## License

[MIT](https://choosealicense.com/licenses/mit/) license, Copyright (c) 2022 Zavolan Lab, Biozentrum, University of Basel 


## Contributers
Eric Boittier, Bastian Wagner, Quentin Badolle

## More info:
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






