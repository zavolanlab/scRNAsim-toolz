# Read Sequencer

## Usage
```
usage: read-sequencer [-h] [-i INPUT] [-r READ_LENGTH] [-n N_RANDOM] [-s CHUNK_SIZE] [-v] output

Simulates sequencing of DNA sequences specified by an FASTA file.

positional arguments:
  output                path to FASTA file

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to FASTA file
  -r READ_LENGTH, --read-length READ_LENGTH
                        read length for sequencing
  -n N_RANDOM, --n_random N_RANDOM
                        n random sequences. Just used if inputfasta file is not specified.
  -s CHUNK_SIZE, --chunk-size CHUNK_SIZE
                        chunk_size for batch processing
```

Examples:
```
read-sequencer -i tests/read_sequencer/files/50_seqs_50_1000_bp.fasta -r 100 sequenced_reads.fa
read-sequencer -n 50 -r 100 random_reads.fa
```

## Overview

Read Sequencer is a python package to simulate sequencing. 
It reads fasta files, simulate sequencing with specified read length and writes the resulting sequences into a new fasta file.


## Contributors and Contact Information

Christoph Harmel - christoph.harmel@unibas.ch  
Michael Sandholzer - michael.sandholzer@unibas.ch  
Clara Serger - c.serger@unibas.ch  

