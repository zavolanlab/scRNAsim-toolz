# Structure generator

## Usage

```
usage: structure-generator [-h] [-p PROB_INCLUSION] [--log LOG] [-v] transcripts annotation

positional arguments:
  transcripts           Path to csv file with number of transcripts (ID,Count).
  annotation            Path to gtf-file with exon annotation.

options:
  -h, --help            show this help message and exit
  -p PROB_INCLUSION, --prob-inclusion PROB_INCLUSION
                        Probability of intron inclusion.
  --log LOG             Level of logging. Can be one of ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]
```

Example:

`structure-generator -p 0.3 tests/structure_generator/files/Transcript1.csv tests/structure_generator/files/Annotations2.gtf`
## Overview

The human body contains a countless variety and diversity of cell types, states, and interactions. We wish to understand these tissues and the cell types at much deeper level. Single-cell RNA-seq (scRNA-seq) offers a look into what genes are being expressed at the level of individual cells. Overall this method allows one to identify cell types, find rare or unidentified cell types or states, identify genes that are differently expressed in different cell types, and explore changes in expression whilst including spatial, regulatory, and protein interactions. 

We hope that others would find use for this transcript_structure generator that allows one to take input gtf-files of specific gene transcripts and outputs a gtf-file containing intron/exon structures per input transcript. Moreover, one can specify a probability for intron-inclusion which is used to simulate incorrect splicing. 



Input:
- csv-formatted file ("ID,Count") with counts for individual transcripts
- probability of intron inclusion (float in range [0,1])
- gtf-formatted file with exon coordinates of the transcripts included in the csv file

Output:
- gtf-formatted file containing generated intron/exon structures per transcript
- csv-formatted file ("NewTranscriptID,ID,Count") with
	- id of generated transcript
	- id of original transcript (without intron inclusions)
	- count


To generate the sampled transcripts, run

```
structure-generator --prob-inclusion <probability_inclusion> [--log "INFO"] <transcripts_file> <annotations_file>
```

where the transcripts file should be csv-formatted, the annotation file gtf-formatted and the inclusion probability for introns a float in the range [0,1]. The log parameter is optional and can be one of `["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]`. The default is `INFO`.

Sample Transcripts and Annotation files can be found in the repository under tests/.

# License

MIT license, Copyright (c) 2021 Zavolan Lab, Biozentrum, University of Basel

# Contributers

Larissa Glass  
Michael Zimmermann  
Andri Fraenkl


