# Read Sequencer

## Overview

Read Sequencer is a python package to simulate sequencing. 
It reads fasta files, simulate sequencing with specified read length and writes the resulting sequences into a new fasta file.


## Installation from github 

Read Sequencer requires Python 3.9 or later.

Install Read Sequencer from Github using:

```
git clone https://git.scicore.unibas.ch/zavolan_group/tools/read-sequencer.git
cd read-sequencer
pip install . 
```

## Usage

```
usage: read_sequencer [-h] [-i INPUT] [-r READ_LENGTH] [-n N_RANDOM] [-s CHUNK_SIZE] output 
Simulates sequencing of DNA sequences specified by an FASTA file.

positional arguments:
  output                path to FASTA file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to FASTA file
  -r READ_LENGTH, --read-length READ_LENGTH
                        read length for sequencing
  -n N_RANDOM, --n_random N_RANDOM
                        n random sequences. Just used if input fasta file is not specified.
  -s CHUNK_SIZE, --chunk-size CHUNK_SIZE
                        chunk_size for batch processing

```

## Docker

The docker image is available on docker hub: https://hub.docker.com/r/grrchrr/readsequencer

```
docker pull grrchrr/readsequencer
docker run readsequencer --help
```

## Contributors and Contact Information

Christoph Harmel - christoph.harmel@unibas.ch  
Michael Sandholzer - michael.sandholzer@unibas.ch  
Clara Serger - c.serger@unibas.ch  

