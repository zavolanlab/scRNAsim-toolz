import pandas as pd

def exon_extraction_from_gtf(gtf_filename,output_filename):
    gtf = pd.read_table(gtf_filename,skiprows=5,header=None)
    exons = gtf[gtf[2]=="exon"]
    features = list(exons[8])

    transcript_id_list = []
    gene_id_list = []
    for x in range(len(features)):
        newlist = features[x].split(";")
        transcript_id_list.append(str(newlist[2])[16:-1])
        gene_id_list.append(str(newlist[0])[9:-1])

    bed = {"chr":exons[0],"start":exons[3],"end":exons[4],"transcript_id":transcript_id_list,"score":exons[5],"strand":exons[6],"gene_id":gene_id_list}
    bed = pd.DataFrame(bed)
    bed.to_csv(output_filename,sep="\t",index=False)

<<<<<<< HEAD

bed = {"chr":exons[0],"start":exons[3],"end":exons[4],"transcript_id":superlist,"score":exons[5],"strand":exons[6],"gene_id":idlist}
class bed:
   def__init__(self, exons, chr, start, end, transcript_id, score, strand, gene_id):
       self.exons = exons
       self.chr = exons[0]
       self.start = exons[3]
       self.end = exons[4]
       self.transcript_id = superlist
       self.score = exons[5]
       self.strand = exons[6]
       self.gene_id = idList

"""Creates BED from GTF for bedtools.

    This class defines a BED from exon annotation from a GTF, for use in bedtools to get sequences with transcript ID as header.
    Parameters
    ----------
    arg1 : GTF file.

    Returns
    -------
    Class
        A class which defines columns in standard BED format.



    Raises
    ------
    TypeError
        ValueError: Not all columns found in GTF.
    """
bed = pd.DataFrame(bed)
bed.to_csv("bed_file.bed",sep="\t",index=False)
bed[(bed["gene_id"]=="ENSG00000160072")|(bed["gene_id"]== "ENSG00000142611")|(bed["gene_id"]=="ENSG00000232596")].to_csv("test.bed",sep="\t",index=False,header=None)
=======
    # This line is used to generate a test file from some of the manually selected gene_ids for now (Plans to make it choose randonly in future)
    bed[(bed["gene_id"]=="ENSG00000160072")|(bed["gene_id"]== "ENSG00000142611")|(bed["gene_id"]=="ENSG00000232596")].to_csv("test.bed",sep="\t",index=False,header=None)
>>>>>>> f35b4d0acdb49bce8ccda3f322fb4b5486737b75
