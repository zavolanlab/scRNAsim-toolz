#!/usr/bin/env nextflow

nextflow.enable.dsl=2

```c
/*   
 * Define the input parameters"""Takes as input FASTA file
                                     of cDNA sequences, a CSV/TSV with sequence
                                     counts, and mean and std. dev. of fragment
                                     lengths and 4 nucleotide probabilities
                                     for the cuts. Outputs most terminal
                                     fragment (within desired length range)
                                     for each sequence."""
 */ 
params.fasta_file = "$projectDir/tests/test_files/test,fasta"
params.counts_file = "$projectDir/tests/test_files/test.csv"        
params.sep = "$projectDir/data/yeast/sep/sep.csv"
params.outdir = "results" 


/* Log some information for the user */

log.info """\
         R N A S E Q - N F   P I P E L I N E 
         ===================================       
         fasta_file   : ${params.fasta_file}      
         counts_file  : ${params.counts_file}             
         outdir       : ${params.outdir}          
         """                                        
         .stripIndent() 
         
/* 
 * Define the `file_validation` process:
 * Validate input files exist and are the correct format
 */ 

process file_validation { 

    input:
    path fasta_file
    path counts_file
    path sep 
    
    output:
    tuple: fasta dict and sequence for counts file
    
    script: 
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """ 
}    
    

/*
 * Define the get_cut_number process:
 * Get the number of cuts for a particular sequence 
 */

process get_cut_number {

    tag "get_cut_numbern on $n_cuts" 
    publishDir "${params.outdir}/get_cut_number", mode:'copy' 
    
    input:
    path index
    tuple val(n_cuts), path(seq_len, mean)                                                                                   
    
    output:
    path(n_cuts)
    
    script:
    """
      
    """                                                      
}


/* 
 * Define the fragmentation process:
 * Fragment cDNA sequences and select terminal fragment 
 */

process fragmentation {

    tag "fragmentation on $fasta, seq_counts, nuc_probs, mu_length, std"
    
    input:
    dict val(fasta), pd.DataFrame(seq_counts), dict(nuc_probs),int(mu_length),int(std)
    
    output:
    path("term_frags")
    
    script:
    """
    mkdir fastqc_${sample_id}_logs
    
    """
}



/* Start the job:
 * initialize variables
 */
 
Channel
    .fromFilePairs( params.reads, checkIfExists:true ) 
    .set { read_pairs_ch }     


/* The "main" function:
 * Use CLI arguments to fragment sequences and output text file with selected terminal fragments
 */

workflow {
  file_validation_ch=file_validation(params.fasta_file, params.counts_file, params.sep)
  get_cut_number_ch = get_cut_number(seq_len, mean)
  framentation_ch = fregmentation(fasta, seq_counts, nuc_probs, mu_length, std)                                                                           }                                                                                                                     


/* Book keeping upon workflow completion */
workflow.onComplete {
   log.info (workflow.success ? "\nDone! Open the following report in your browser -->     $params.outdir/multiqc/multiqc_report.html\n" : "Oops .. something went wrong") 
   )
}                                                                                                                     
```
