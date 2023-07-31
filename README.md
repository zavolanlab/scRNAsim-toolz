# Transcript Sampler

## Overview
This workflow samples representative transcripts per gene, in proportion to their relative abundance levels. Sampling is done by Poisson sampling. 

This workflow takes as input:
- Path to genome annotation file in gtf format
- Path to csv or tsv file with transcript IDs and expression levels
- Path to output sample gtf file 
- Path to output sample transcript IDs and counts
- Integer of number of transcripts to sample
 
The outputs are :
- trancript sample gtf file
- csv file containing sample transcript IDs and counts.
 
## Installation from github
Transcript sampler requires Python 3.9 or later.

Install Transcript sampler from Github using:

```
git clone https://git.scicore.unibas.ch/zavolan_group/tools/transcript-sampler.git
cd transcript-sampler
pip install . 
```

## Usage
```
usage: transcript-sampler [-h] --input_gtf INPUT_GTF --input_csv INPUT_CSV --output_gtf OUTPUT_GTF --output_csv OUTPUT_CSV --n_to_sample N_TO_SAMPLE

Transcript sampler

options:
  -h, --help            show this help message and exit
  --input_gtf INPUT_GTF
                        GTF file with genome annotation (default: None)
  --input_csv INPUT_CSV
                        CSV or TSV file with transcripts and their expression level (default: None)
  --output_gtf OUTPUT_GTF
                        Output path for the new GTF file of representative transcripts (default: None)
  --output_csv OUTPUT_CSV
                        Output path for the new CSV file of representative transcripts and their sampled number (default: None)
  --n_to_sample N_TO_SAMPLE
                        Total number of transcripts to sample (default: None)
```

Example : 

`transcript_sampler --input_gtf="tests/input_files/test.gtf" --input_csv="tests/input_files/expression.csv" --output_gtf="output_files/output.gtf" --output_csv="output_files/output.csv" --n_to_sample=100`
