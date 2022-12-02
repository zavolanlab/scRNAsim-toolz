# Transcript Sampler

This workflow takes as input:
 - genome annotation gtf file
 - expression levels of each gene
 - csv file with transcript IDs and expression levels
 
 The output is a trancript sample gtf file and csv file containing transcript IDs and counts.
 
 The workflow can be run via the command line as
 
 python scripts/exe.py --annotation {gtf input file} --output_csv {output csv file} --transcript_number {number of transcripts} --output_gtf {output gtf file} --input_csv {input csv file}
