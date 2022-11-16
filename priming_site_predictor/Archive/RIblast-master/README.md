# RIblast
RIblast is ultrafast RNA-RNA interaction prediction software based on seed-and-extension algorithm for comprehensive lncRNA interactome analysis.

## Version
Version 1.2.0 (2019/11/02)

## Acknowledgements
We used BL* energy model as RNA secondary structure energy model.
You can download this energy model from  http://www.cs.ubc.ca/labs/beta/Projects/RNA-Params
We thank Dr. Mirela Andronescu for the development of this model.

## Usage
RIblast consists of database construction search step and rna interactin search step. Firstly, you have to generate formatted database files from FASTA fomatted RNA sequences file using RIblast "db" command. RIblast "db" command requires 2 options ([-i InputFastaFile] and [-o OutputDbName]). Then, you can search RNA-RNA interaction of a query seuence to database sequences using RIblast "ris" command. RIblast "ris" command requires 3 options ([-i InputFastaFile], [-o OutputFileName] and [-d DatabaseFileName]). 

## Example
    ./RIblast db -i dbRNA.fa -o test_db
    ./RIblast ris -i queryRNA.fa -o output.txt -d test_db

## Command and Options
    db: convert a FASTA file to RIblast database files  

     RIblast db [-i InputFastaFile] [-o OutputDbName] [-r RepeatMaskingStyle]  
                [-s LookupTableSize] [-w MaximalSpan] [-d MinAccessibleLength]  
   
    Options:
    (Required)
        -i STR    RNA sequences in FASTA format
        -o STR    The database name
        
    (Optional) 
        -r INT    Designation of repeat masking style 0:hard-masking, 1:soft-masking, 2:no-masking [default:0]
        -s INT    Lookup table size of short string search [default: 8]
        -w INT    The constraint of maximal distance between the bases that form base pairs. This parameter have to be set to 20 and over. [default: 70]
        -d INT    Minimum accessible length for accessibility approximation [defualt:5]
        
    ris: search RNA-RNA interaction between a query and database sequences
    
    RIblast ris [-i InputFastaFile] [-o OutputFileName] [-d DatabaseFileName]
                [-l MaxSeedLength] [-e HybridizationEnergyThreshold] [-f InteractionEnergyThreshold]
                [-x DropOutLengthInGappedExtension] [-y DropOutLengthInUngappedExtension]
                [-g OutputEnergyThreshold] [-s OutputStyle]
                
    Options:
    (Required)
        -i STR    RNA sequences in FASTA format
        -o STR    Output file name
        -d STR    The database name
        
    (Optional)
        -l INT    Max size of seed length [default:20]
        -e DBL    Hybridization energy threshold for seed search [default: -6.0]
        -f DBL    Interaction energy threshold for removal of the interaction candidate before gapped extension [default: -4.0]
        -x INT    DropOut Length in gapped extension [defualt:16]
        -y INT    DropOut Length in ungapped extension [defualt:5]
        -g DBL    Energy threshold for output [defualt:-8.0]
        -s INT    Designation of output format style 0:simplified output style, 1:detailed output style [defualt:0]

## Output file format
RIblast outputs detected interactions as follows.
An interaction is expressed in five columns.The first, second, third, fourth and fifth column of an interaction describes intearction id, name of query RNA, name of target RNA, interaction energy of the interaction, and interacted regions ([region in query]:[region in target]), respectively.

Example:  

RIblast ris result  
input:test_input.txt,database:test_db,RepeatFlag:0,MaximalSpan:100,MinAccessibleLength:3,MaxSeedLength:20,InteractionEnergyThreshold:0,HybridEnergyThreshold:3,FinalThreshold:0,DropOutLengthWoGap:5,DropOutLengthWGap:16  
Id,Query name, Query Length, Target name, Target Length, Accessibility Energy, Hybridization Energy, Interaction Energy, BasePair  
0,qrna,100,target_rna1,200,9.44198,-19.66,-10.218,(4-21:30-13)  
1,qrna,100,target_rna2,150,5.95465,-15.06,-9.10535,(72-83:185-170)  

Query Length and Target Length means the region of non-repeat region when repeat sequences are masked.  
If you need the information of interacted base pair, please designate the output format style as detailed output style.

## License
This software is released under the MIT License, see LICENSE.txt.

## Changelogs
2019/11/02 Version 1.2.0 was released. The calculation speed is about twice as fast as before.  
2017/11/30 Version 1.1.3 buf fix:calculation of dangling energy.  
2017/07/06 Version 1.1.2 I changed output file format.  
2016/12/07 Version 1.1.1 I changed output file format.  
2016/11/21 Version 1.1.0 was released.  
2016/10/10 Version 1.0.2 bug fix: I fixed a bug in accessibility calculation.  
2016/09/27 Version 1.0.1 bug fix: I fixed a bug in loading fastafile and calculation of dangling energy. I also changed output file format.  
2016/08/31 Version 1.0.0 was released.

## Reference
Tsukasa Fukunaga and Michiaki Hamada. "RIblast: An ultrafast RNA-RNA interaction prediction system based on a seed-and-extension approach." btx287, Bioinformatics (2017)
