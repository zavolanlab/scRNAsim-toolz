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
