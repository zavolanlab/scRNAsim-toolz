# Terminal fragment selector

## Usage
```
usage: fragment-selector [-h] --fasta FASTA --counts COUNTS -o OUTPUT [--mean MEAN] [--std STD] [-s SIZE] [--sep SEP] [-v]

Takes as input FASTA file of cDNA sequences, a CSV/TSV with sequence counts, and mean and std. dev. of fragment lengths and 4 nucleotide probabilities for the cuts. Outputs most terminal fragment (within desired length
range) for each sequence.

options:
  -h, --help            show this help message and exit
  --fasta FASTA         Path to FASTA file with cDNA sequences
  --counts COUNTS       Path to CSV/TSV file with sequence counts
  -o OUTPUT, --output OUTPUT
                        output file path
  --mean MEAN           Mean fragment length (default: 300)
  --std STD             Standard deviation fragment length (defafult: 60)
  -s SIZE, --size SIZE  Chunk size for batch processing
  --sep SEP             Sequence counts file separator.
```

Example:

`fragment-selector --fasta tests/fragment_selector/files/test.fasta --counts tests/fragment_selector/files/test.csv --mean 50 --output fragments.fa`
## Overview
Simulating single cell RNA library generation (scRNA-seq)

This repository is as part of the Uni Basel course <E3: Programming for Life Science – 43513>. To test the accuracy of scRNA-seq data we generated the *synthetic data*. That is, we reconstruct the properties of the experimental data set and determine whether the computational analysis can recover properties of the data that was assumed in the simulation. This is never trivial since setting the ground truth is much needed in the computational method to evaluate the result. 


As part of the sub-project, we implemented python code for selecting terminal fragments. Detailed distribution used for the selecting fragments can be found below, summarised in [this paper](https://www.nature.com/articles/srep04532#MOESM1).
> Next Generation Sequencing (NGS) technology is based on cutting DNA into small fragments and their massive parallel sequencing. The multiple overlapping segments termed “reads” are assembled into a contiguous sequence. To reduce sequencing errors, every genome region should be sequenced several dozen times. This sequencing approach is based on the assumption that genomic DNA breaks are random and sequence-independent. However, previously we showed that for the sonicated restriction DNA fragments the rates of double-stranded breaks depend on the nucleotide sequence. In this work we analyzed genomic reads from NGS data and discovered that fragmentation methods based on the action of the hydrodynamic forces on DNA, produce similar bias. Consideration of this non-random DNA fragmentation may allow one to unravel what factors and to what extent influence the non-uniform coverage of various genomic regions.

As a whole project, we implemented a procedure for sampling reads from mRNA sequences, incorporating a few sources of “noise”. These include the presence of multiple transcript isoforms from a given gene, some that are incompletely spliced, stochastic binding of primers to RNA fragments and stochastic sampling of DNA fragments for sequencing. We will then use standard methods to estimate gene expression from the simulated data. We will repeat the process multiple times, each time corresponding to a single cell. We will then compare the estimates obtained from the simulated cells with the gene expression values assumed in the simulation. We will also try to explore which steps in the sample preparation have the largest impact on the accuracy of gene expression estimates.



CLI arguments:
- fasta (required): Path to FASTA file with cDNA sequences
- counts (required): Path to CSV/TSV file with sequence counts
- output (required): Output file path

- mean: Mean fragment length (default: 300)
- std: Standard deviation fragment length (default: 60)
- size: Chunk size for batch processing (default: 10000)
- sep: Sequence counts file separator (default: ",")

"""Takes as input FASTA file of cDNA sequences, a CSV/TSV with sequence counts, and mean and std. dev. of fragment lengths. Outputs most terminal fragment (within desired length range) for each sequence."""

Output:
- Text file with most terminal fragments for each sequence.


# License

MIT license, Copyright (c) 2021 Zavolan Lab, Biozentrum, University of Basel

# Contact
zavolab-biozentrum@unibas.ch

