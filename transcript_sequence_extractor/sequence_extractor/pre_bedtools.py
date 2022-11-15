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

    # This line is used to generate a test file from some of the manually selected gene_ids for now (Plans to make it choose randonly in future)
    bed[(bed["gene_id"]=="ENSG00000160072")|(bed["gene_id"]== "ENSG00000142611")|(bed["gene_id"]=="ENSG00000232596")].to_csv("test.bed",sep="\t",index=False,header=None)
