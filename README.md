# Transcript Sampler

This workflow samples representative transcripts per gene, in proportion to their relative abundance levels. Sampling is done by Poisson sampling. 

**This workflow takes as input:**
 - Path to genome annotation file in gtf format
 - Integer of number of transcripts to sample
 - Path to csv or tsv file with transcript IDs and expression levels
 - Path to output sample gtf file 
 - Path to output sample transcript IDs and counts
 
 **The outputs are :**
 - trancript sample gtf file
 - csv file containing sample transcript IDs and counts.
 
 **The workflow can be run via the command line as**
 
 `python transcript_sampler/new_exe.py --input_gtf={gtf input file} --input_csv={input csv file} --output_gtf={output gtf file} --output_csv={output csv file} --n_to_sample={number of transcripts}`

 Example : 

 `python transcript_sampler/new_exe.py --input_gtf="input_files/test.gtf" --input_csv="input_files/expression.csv" --output_gtf="output_files/output.gtf" --output_csv="output_files/output.csv" --n_to_sample=100`
