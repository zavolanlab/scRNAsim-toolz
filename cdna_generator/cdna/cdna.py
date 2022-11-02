import sys

def translate(res):
    translate_dict = {"A": "T", "U": "A", "G": "C", "C":"G"}
    if res not in translate_dict.keys():
        print("cDNA residue not A,T,U or G ")
        sys.exit(1)
    return translate_dict[res]

class cDNA_Gen:
    def __init__(self, 
                 fasta, 
                 gtf, 
                 cpn, 
                 output_fasta = "cDNA.fasta",
                 output_csv = "cDNA.csv"
                ):
        # inputs
        self.output_fasta = output_fasta
        self.output_csv = output_csv
        
        # variables
        self.prime_sites = []
        self.fasta_seq = ""
        
        self.main()
        
    def main(self):
        self.read_fasta()
        self.read_gtf()
        
    def order_priming_sites(self):
        pass
    
    def generate_cdna(self):
        pass
    
    def read_fasta(self):
        pass
    
    def read_gtf(self):
        pass
    
    def write_fasta(self):
        pass
    
    def read_copy_numbers(self):
        pass
    
    def return_output(self):
        return self.output_fasta, self.output_csv
    

    
    
    
