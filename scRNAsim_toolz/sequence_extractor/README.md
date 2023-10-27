# Sequence Extractor

## Usage
```
usage: sequence-extractor [-h] --mode {pre_bedtools,post_bedtools} [-i INPUT_FASTA_FILE] [-o OUTPUT_FILE_NAME] [-p POLY_A_LENGTH] [--input-gtf-file INPUT_GTF_FILE] [--output-bed-file OUTPUT_BED_FILE] [-v]

extracts transcript sequences from genome sequence andouputs transcripts with PolyA tail added to them

options:
  -h, --help            show this help message and exit
  --mode {pre_bedtools,post_bedtools}
                        Select the mode of operation('pre_bedtools' or 'post_bedtools').
  -i INPUT_FASTA_FILE, --input-fasta-file INPUT_FASTA_FILE
                        Fasta-formatted file obtained from bedtools
  -o OUTPUT_FILE_NAME, --output-file-name OUTPUT_FILE_NAME
                        Name of the output fasta file
  -p POLY_A_LENGTH, --polyA-length POLY_A_LENGTH
                        Length of the polyA tail to be added (def: 250)
  --input-gtf-file INPUT_GTF_FILE
                        Ordered and processed gtf file for 'pre_bedtools' mode.
  --output-bed-file OUTPUT_BED_FILE
                        Bed file with only exons with strandednesstaken into account for 'pre_bedtools' mode.
```

Example:
```
sequence-extractor --mode pre_bedtools --input-gtf-file tests/sequence_extractor/files/test.gtf --output-bed-file tests/sequence_extractor/files/output.bed
gunzip -k tests/sequence_extractor/files/human.chr1.fa.gz 
bedtools getfasta -fi tests/sequence_extractor/files/human.chr1.fa -bed tests/sequence_extractor/files/output.bed -name -s -fo tests/sequence_extractor/files/output.fasta
sequence-extractor --mode post_bedtools --input-fasta-file tests/sequence_extractor/files/output.fasta --polyA-length 250 --output-file-name polyA_output.fa

sequence-extractor --mode post_bedtools --input-fasta-file tests/sequence_extractor/files/post_bedtools_test.fa --polyA-length 250 --output-file-name polyA_output.fa
```

## Overview 

Given a gtf specification of transcript exon/intron structures and the genome sequence, construct the nucleotide sequence of the transcripts and add poly(A) tails.

__Input:__

* Gtf file with exon/intron structures of transcripts
* File with genome sequence
* Length of the poly(A) tail
* Dictionary of expected nucleotide frequencies in poly(A) tail


__Output:__

For each transcript, the list of exons should be traversed from 5' to 3', the sequences of the exons need to be extracted from the genome given the coordinates and then pasted together. At the end, a tail of the specified length should be added at the 3' end of the transcript, given a vector of mono-nucleotide frequencies (of course, the frequency of A's will be much higher than of any other nucleotide).


## Design plan

### 1-  Obtain gtf file and also generate a test file for code validation (sampled transcript gtf, from Group 2) :  
a. To use as main test file: Reference gtf file  
b. Subset the gtf file so that only the rows labelled exon remain  
c. Extract gtf file data to make a BED file (bedtools getfasta -split only works on BED12 files)  
d. For the sample test file, use just one gene based on gene id. Later on, the sample test file can be expanded as required

### 2- To run bedtools on the genome with data from sample/final bed file and extract exon sequences for all transcripts (bedtools getfasta) :
a. To use as genome sequence: https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz (Human Genome release 107 from Ensembl)  
b. Using bedtools getfasta and its options, get the output
c. Output will be: transcript sequences containing all exon data in 5’ to 3’ direction without poly A tail  

### 3-Function to add poly A tail at the 3’ end (right side) of RNA sequence:  
a. Input: The output transcripts from the previous step
b. This is based on length of tail (sample a length of 250 nucleotides)  
c. It needs to take into account the probability (weight/relative frequency) of each nucleotide (with A having the largest frequency defined)- this will be saved as a dictionary. {‘A’:freqA,‘U’:freqU,‘G’:freqG,‘C’:freqC}  
d. Append via e.g. str.join(), str.ljust()  
e. Output the final transcript sequences as a .fasta file. (Final Output)  


## License

[MIT](https://choosealicense.com/licenses/mit/) license, Copyright (c) 2022 Zavolan Lab, Biozentrum, University of Basel 


## Contributers
Samuel Mondal, Ahmed Hassan Hussein H.Mahmoud, Gina Boot

