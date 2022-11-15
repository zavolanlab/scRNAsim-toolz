import pandas as pd

gtf = pd.read_table("Homo_sapiens.GRCh38.107.gtf.gz",skiprows=5,header=None)


exons = gtf[gtf[2]=="exon"]
feat = list(exons[8])
superlist = []
idlist = []
for x in range(len(feat)):
    newlist = feat[x].split(";")
    superlist.append(str(newlist[2])[16:-1])
    idlist.append(str(newlist[0])[9:-1])


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
