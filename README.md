# scRNAsim: Simulating single cell RNA (scRNA-seq) library generation 
The projects implements a simulation of single cell RNA sequencing (scRNA-seq), accounting for some common sources noise that complicate the analysis of the resulting data.

### Setting up the virtual environment

Create and activate the environment with necessary dependencies with Conda:

```bash
conda env create -f environment.yml
conda activate scrnasim-toolz
```

### Tools

The tools available in this repo are:
1. Transcript sampler
2. Structure generator
3. Sequence extractor
4. Priming site predictor
5. cDNA generator
6. Fragment selector
7. Read sequencer

### Description

Although all cells in a multicellular organism carry the same genomic information, they differ a lot in their function, due to the fact that they are equipped with distinct toolboxes of molecular functions, implemented by different proteins and RNAs. Thus, being able to detect and measure the abundance of gene products (RNAs and/or proteins) in individual cells holds the key to understanding how organisms are organized and function. In the past decade, much progress has been made in the development of technologies for single cell RNA sequencing. They make use of microfluidic devices that allow RNA-seq sample preparation for individual cells encapsulated in droplets, followed by pooling of the resulting DNA fragments and sequencing. The broadly used 10x Genomics technology uses oligo-dT primers to initiate the cDNA sequencing from the poly(A) tails of fragmented RNAs. Subsequent sequencing yields *libraries* of relatively short (100-200 nucleotides) *reads* that come predominantly from the 3’ ends of RNAs given the priming on the poly(A) tail. As in the ideal case (no amplification bias) each read came from the end of one mRNA, simply counting the reads that map to mRNAs of individual genes provides estimates of the expression levels of those genes within the respective cell. Currently, typical data sets cover thousands of genes in tens-to-hundreds of thousands of cells. However, we are still far from being able to prepare ideal libraries, for many reasons. First, as gene expression is a bursty, stochastic process, there will be fluctuations in the number of RNAs (corresponding to a given gene) that are present in any one cell at the time of sampling, even when the time-average of those RNA numbers were to be the same across cells. Secondly, the sample preparation steps are carried out by various enzymes with limited efficiency. This leads to substantial fluctuations in the number of molecules that are “captured” for a gene in a given cell, even if all cells were to have the same abundance of these molecules at the time of sampling. Third, the biochemical reactions that are part of sample preparation do not have absolute specificity. A clear example is the priming of the cDNA synthesis with oligo(dT): although the primer is intended for the poly(A) tails at the 3’ ends of RNAs, it is clear that the primer also binds to A-rich stretches that are internal to transcripts, and especially located in intronic regions. Finally, a conceptual issue with the single cell data is that one cannot apply the principle of averaging measurement values across replicate experiments to obtain more precise estimates, because we do not know which cells could be considered replicates of each other (if that is at all conceivable).

For all of these reasons, testing the accuracy of computational analysis methods for scRNA-seq data is not trivial. One generally does not have “ground truth” data on which to benchmark computational methods, but there are various ways in which scientists get around this problem. One is to generate *in vitro* data sets for which the analysis should produce results that fall within predictable bounds. For example, this was done in one case as follows (https://www.nature.com/articles/nmeth.2930): after pooling the RNAs from many cells, the authors have sampled “single-cell equivalent” pools of RNAs, and passed them through the sample preparation protocol. In this case, the assumption that the single-cell equivalents are noisy samples from a unique vector of gene expression levels holds. Fluctuations in the estimates of expression levels were therefore due to 2 factors only: the sampling of the single-cell equivalents and the sampling steps in the sample preparation protocol.
The second approach to assessing the accuracy of computational analysis methods is to use *synthetic data*. That is, to generate data sets by simulating the experimental steps and determine whether the computational analysis can recover properties of the data that was assumed in the simulation. For example, if the goal of the computational analysis is to infer gene expression levels from scRNA-seq data, then one simulates such data assuming specific transcript abundances, which should be recovered by the computational method. In general, it is very difficult to accurately model each step of the experimental procedure, and therefore, simulations still leave out some (possibly a lot) of the complexities of the experiment. Thus, the fact that a computational method performs well on simulated data provides more of a sanity check on the method than the confidence that the method will give accurate results on *real* data. Nevertheless, such sanity checks should be done. Furthermore, simulations can help build intuitions as to which steps of the experiment have the largest consequences for the outcome, where specific behaviors may come from etc.

In this project we will implement a procedure for sampling reads from mRNA sequences, incorporating a few sources of “noise”. These include the presence of multiple transcript isoforms from a given gene, some that are incompletely spliced, stochastic binding of primers to RNA fragments and stochastic sampling of DNA fragments for sequencing. We will then use standard methods to estimate gene expression from the simulated data. We will repeat the process multiple times, each time corresponding to a single cell. We will then compare the estimates obtained from the simulated cells with the gene expression values assumed in the simulation. We will also try to explore which steps in the sample preparation have the largest impact on the accuracy of gene expression estimates.

Inputs to the simulation:
1. Total number of transcripts per cell and csv-formatted table “GeneID,Counts” specifying the number of transcripts that are expressed, on average, for gene GeneID in a given cell type. These can come for example from a bulk RNA-seq experiment of sorted cells of a given type. 
2. File with the genome sequence
3. gff/gtf-formatted file with the transcript annotation of the genome
4. Output directory
5. Number of reads to sequence
6. Number of cells to simulate
7. Mean and standard deviation of RNA fragment length
8. Read length and nucleotide frequencies for padding
9. Probability of intron inclusion - considered constant per intron to start with, can be extended to intron-specific. In the latter case, estimates could be obtained from bulk RNA-seq data by dividing the average per-position coverage in a given intron by the average per-position coverage of the gene, or of flanking exons.
10. Option to add poly(A) tails to transcripts and an associated function for generating these tails (with specific length distribution and non-A nucleotide frequency).
11. Parameters for evaluating internal priming: primer sequence, function implementing the constraints on priming sites (accessibility, energy of interaction, perfect matching at last primer position etc.).

The workflow will then go through the following steps, repeated for the *number of cells* times, to simulate data for multiple cells.
1. Pick the number of transcripts coming from each gene. As #input 1 we get a file with the expression level of individual transcripts from some real sample. For simplicity, we first pick a representative transcript per gene, e.g. with most annotation support (support level 1 or TSL=1). Then, given a total number of transcripts per cell (input #1, we generate, for each representative transcript, a Poisson sample given the average count from input #1.
2. Generate the exon/intron structure for each transcript. For each gene there is a reference set of exons (specified in input #3), but to account for the possibility that the transcript is not completely processed, introns are included in individual transcripts. This is done by going through all the possible introns of a transcript and choosing which ones to include (according to input #9). Then the generated structures are written to a gtf file (because new exons are effectively generated which look like exon_n;intron_n;exon_n+1) and the number of transcripts with each unique structure is also saved.
3. To start the "sample preparation" process, the transcript sequences need to be constructe based on the structure specified in the gff/gtf file and the genome sequence. Poly(A) tails are added to transcripts accoroding to input #10.
4. Knowing that during the cDNA synthesis process poly(A) stretches act as priming sites, we predict how likely it is to initiate synthesis at every position on each transcript. We assume that this depends on the energy of binding between the poly(T) primer and the stretch of transcript starting at the specified position. We run an external program to predict this energy and we thus obtain, for each each transcript a list of positions where priming as non-0 probability (according to hybridization parameters specified by input #11).
5. The possible priming sites are sampled with the probabilities computed at the previous step, to pick a site for generating the complementary DNA.
6. The resulting cDNAs are fragmented according to the parameters specified by input #7, and the end-fragments of the cDNA are selected.
7. The terminal fragments from the previous step are sampled according to input #5, to pick a fragment for sequencing. Then a piece of length input #8 is taken fromm the 5' end of the fragment to form a read. If the fragment is shorter than the read length (input #8), the fragment is padded with random sequence, given a vector of relative probability for A,C,G,T to appear in the random sequence (input #8). The output of this step will be a fasta file with "sequenced reads", which is the output of the simulation.

We can use the simulated data to evaluate methods that map reads to genome/transcripts and quantify transcript abundance. For instance, the following analyses can be implemeted.
1. Map reads to genome - use the STAR aligner to map the reads obtained from each cell to the genome, obtaining one bam file per cell.
2. Quantify gene expression in every cell - assign reads to genes based on the genome anotation (taking into account that reads can map to more than one location in the genome) and generate a table of GeneID-count_in_cell_1-count_in_cell_2-count_in_cell_3-...
3. Calculate mean and variance of expression level for each gene across cells
4. Visualize these numbers as "M-A plots", which are countour plots where the x-axis is the average expression level, the y-axis is the variance and the contours show the density of genes with specific mean-variance values.
5. Given a GeneID and the counts per cells, construct the histogram of the counts for the gene across cells
6. Check the accuracy of gene expression inference by comparing the average input expression (input #1) with the inferred gene expression levels (calculated at step #3 above).

Each processing step should be wrapped into a *Nextflow* process, and the whole simulation should be put together as a *Nextflow* workflow.

For each step, test cases should be generated.
