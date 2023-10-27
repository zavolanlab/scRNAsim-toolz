# Priming Site Predictor of Transcript Sequences

## Usage
```
usage: priming-site-predictor [-h] [-f FASTA_FILE] [-p PRIMER_SEQUENCE] [-e ENERGY_CUTOFF] [-r RIBLAST_OUTPUT] [-o OUTPUT_FILENAME] [-v]

Compute potential priming sites using RIBlast.

options:
  -h, --help            show this help message and exit
  -f FASTA_FILE, --fasta-file FASTA_FILE
                        Fasta-formatted file of transcript sequences
  -p PRIMER_SEQUENCE, --primer-sequence PRIMER_SEQUENCE
                        Primer sequence
  -e ENERGY_CUTOFF, --energy-cutoff ENERGY_CUTOFF
                        Energy cutoff for interactions
  -r RIBLAST_OUTPUT, --riblast-output RIBLAST_OUTPUT
                        Path to RIBlast output file
  -o OUTPUT_FILENAME, --output-filename OUTPUT_FILENAME
                        Path where the output gtf should be written
```

Example:
```
priming-site-predictor --riblast-output tests/priming_site_predictor/files/RIBlast_output_example.txt --output-filename priming_sites.gtf
```

## Overview
Priming Site Predictor which uses a seed-and-extension algorithm (*RIblast*: https://github.com/fukunagatsu/RIblast) to *Predict Priming Sites* of oligo dT primers in target sequences. Furthermore, *Binding Energies* are calculated and classified with a threshold value. Additionally, the binding sites are associated with *Binding Probabilities* and stored in a *gtf file* for further processes.

## Usage
```
usage: primingsitepredictor [-h] [-f FASTA_FILE] [-p PRIMER_SEQUENCE] [-e ENERGY_CUTOFF] [-r RIBLAST_OUTPUT] [-o OUTPUT_FILENAME]

Compute potential priming sites using RIBlast.

options:
  -h, --help            show this help message and exit
  -f FASTA_FILE, --fasta-file FASTA_FILE
                        Fasta-formatted file of transcript sequences
  -p PRIMER_SEQUENCE, --primer-sequence PRIMER_SEQUENCE
                        Primer sequence
  -e ENERGY_CUTOFF, --energy-cutoff ENERGY_CUTOFF
                        Energy cutoff for interactions
  -r RIBLAST_OUTPUT, --riblast-output RIBLAST_OUTPUT
                        Path to RIBlast output file
  -o OUTPUT_FILENAME, --output-filename OUTPUT_FILENAME
                        Path where the output gtf should be written
```

## Example
```
primingsitepredictor -f tests/test_files/test_fasta.fasta -p tests/test_files/primer1.fasta -e 0.002 -r tests/test_files/RIBlast_output_example.txt -o tests/test_output.gtf
```

## License
This software is released under the MIT License, see LICENSE.txt.

## Changelogs
2022/11/15 Version 0.1.0 was released.

## Contributors
Max Bär, Sophie Schnider, Robin Christen (University of Basel)

## Acknowledgements
We used the RIblast algorithm created by Tsukasa Fukunaga (https://github.com/fukunagatsu). 

## Reference
Tsukasa Fukunaga and Michiaki Hamada. "RIblast: An ultrafast RNA-RNA interaction prediction system based on a seed-and-extension approach." btx287, Bioinformatics (2017)
