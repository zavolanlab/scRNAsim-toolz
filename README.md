# Transcript Sampler

This workflow sample representative transcripts per gene, in proportion to their relative abundance levels. Sampling is done by poisson sampling. 

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
 
 `python scripts/new-exe.py --annotation {gtf input file} --output_csv {output csv file} --transcript_number {number of transcripts} --output_gtf {output gtf file} --input_csv {input csv file}`

 Exemple : 

 `python scripts\new_exe.py --annotation "input_files\test.gtf" --output_csv "output_files\output_csv.txt" --transcript_number 50  --output_gtf "output_files\output_gtf.gtf" --input_csv "input_files/expression.csv"`


