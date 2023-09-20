# cDNA Generator module
Generate cDNA based on mRNA transcript sequences and the coresponding priming probabilities. 

## Example usage
A simple example can be run from the test_files directory:

    python ../cdna/cli.py -ifa yeast_example.fa -icpn copy_number_input.csv -igt Example_GTF_Input.GTF -ofa cDNA.fasta -ocsv cDNA.csv

## Installation 

    pip install .

## Docker
A docker image is available, to fetch this image:

    docker pull ericdb/my-image

To run a simple example using this image:

    docker run my-image python cdna/cli.py -ifa test_files/yeast_example.fa -icpn test_files/copy_number_input.csv -igt test_files/Example_GTF_Input.GTF -ofa test_files/cDNA.fasta -ocsv test_files/cDNA.csv

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






