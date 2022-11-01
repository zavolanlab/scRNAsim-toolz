# Transcript Sequence Extractor 
### Design plan:
 1. Obtain gtf file and also generate a test file for code validation(sampled transcript gtf, from Group 2) :   
    a. To use as main test file: Reference gtf file  
    b. Subset the gtf file so that only the rows labelled exon remain  
    c. Extract gtf file data to make a BED file (bedtools getfasta -split only works on BED12 files)  
    d. For the sample test file, use just one gene based on gene id. Later on, the sample test file can be expanded as required  
 2. To run bedtools on the genome with data from sample/final bed file and extract exon sequences for all transcripts (bedtools getfasta) :   
    a. To use as genome sequence: https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz (Human Genome release 107 from Ensembl)  
    b. Using bedtools getfasta and its options, get the output  
    c. Output will be: transcript sequences containing all exon data in 5’ to 3’ direction without poly A tail  
 3. Function to add poly A tail at the 3’ end (right side) of RNA sequence:  
    a. Input: The output transcripts from the previous step  
    b. This is based on length of tail (sample a length of 250 nucleotides)  
    c. It needs to take into account the probability (weight/relative frequency) of each nucleotide (with A having the largest frequency defined)- this will be saved as a dictionary. {‘A’:freqA,‘U’:freqU,‘G’:freqG,‘C’:freqC}   
    d. Append via e.g. str.join(), str.ljust()  
    e. Output the final transcript sequences as a .fasta file. (Final Output)
