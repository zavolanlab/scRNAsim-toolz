#!/usr/bin/env nextflow

nextflow.enable.dsl=2

```c
/*   
 * Define the input parameters"""Takes input from csv-formatted file ("ID,Count") 
                                     with counts for individual transcripts. 
                                     Also a gtf-formatted file with exon coordinates 
                                     of the transcripts included in the csv file. 
                                     Outputs gtf-formatted file containing 
                                     generated intron/exon structures per transcript.
                                     Outputs csv-formatted file 
                                     ("NewTranscriptID,ID,Count") with id of 
                                     generated transcript, id of original 
                                     transcript (without intron inclusions) and
                                     Counts. 
 */ 
params.cvs_file = "$projectDir/tests/ #Input cvs file from group 1
params.gtf_file = "$projectDir/tests/ #Input gtf file from group 1       
params.outdir = "results" 


/* Log some information for the user */

log.info """\
         R N A S E Q - N F   P I P E L I N E 
         ===================================       
         cvs_file     : ${params.cvs_file}      
         gtf_file     : ${params.gtf_file}             
         outdir       : ${params.outdir}          
         """                                        
         .stripIndent() 
         
/* 
 * define the `file_validation` process
 * given that the input files exist and are the correct format
 */ 

process file_validation { 

    input:
    path cvs_file
    path gtf_file 
    
    output:
    path cvs_file into struc_gen
    path gtf_file into struc_gen
    
    script: 
    """
    
    """ 
}  

/* 
 * Define the transcript structure generation
 * Read gtf files to define exons/introns and add random mutation coefficient. 
 */

process struc_gen {
    
    input:
    
    
    output:
    
    
    script:
    """
    
    
    """
}



/* Start the job:
 * initialize variables
 */
 
Channel
    .fromFilePairs( params.reads, checkIfExists:true ) 
    .set { read_pairs_ch }     


/* The "main" function:
 * Use CLI arguments to create structure sequences and output text file & gtf file with selected exon structures
 */

workflow {
  file_validation_ch = file_validation(params.cvs_file, params.gtf_file)
  struc_gen_ch = struc_gen( # No idea yet)                                                                           }                                                                                                                     


/* Book keeping upon workflow completion */
workflow.onComplete {
   log.info (workflow.success ? "\nDone! Open the following report in your browser -->     $params.outdir/multiqc/multiqc_report.html\n" : "Oops .. something went wrong") 
   )
}                                                                                                                     
```
