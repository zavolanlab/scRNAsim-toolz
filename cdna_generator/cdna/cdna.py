import sys


def translate(res):
    translate_dict = {"A": "T", "U": "A", "G": "C", "C": "G"}
    if res not in translate_dict.keys():
        print("cDNA residue not A,T,U or G ")
        sys.exit(1)
    return translate_dict[res]


class cDNA_Gen:
    def __init__(
        self, fasta, gtf, cpn, output_fasta="cDNA.fasta", output_csv="cDNA.csv"
    ):
        # inputs
        self.fasta = fasta
        self.gtf = gtf
        self.cpn = cpn
        self.output_fasta = output_fasta
        self.output_csv = output_csv
        # variables
        self.prime_sites = []
        self.fasta_seq = ""
        self.fasta_id = ""
        self.copy_numbers = {}

        self.run()

    def run(self):
        self.read_fasta()
        self.read_gtf()

    def order_priming_sites(self):
        pass

    def generate_cdna(self):
        pass

    def read_fasta(self):
        fasta = open(self.fasta).readlines()
        self.fasta_id = fasta[0]
        print(fasta[0])
        self.fasta_seq = "".join([_.rstrip() for _ in fasta[1:]])

    def read_gtf(self):
        with open(self.gtf) as gtf_file:
            gtf_lines = gtf_file.readlines()
            for line in gtf_lines[:1000]:
                if not line.startswith("#"):
                    temp_gtf = GTF_entry(line)
                    temp_gtf.set_sequence(self.fasta_seq)
                    self.prime_sites.append(temp_gtf)

    def write_fasta(self):
        pass

    def read_copy_numbers(self):
        with open(self.cpn) as cpn_file:
            cpn_lines = cpn_file.readlines()
            for line in cpn_lines:
                csv = line.split(",")
                trans_id = csv[0]
                if trans_id:
                    gene_id = csv[1]
                    count = csv[2]
                    self.copy_numbers[gene_id] = count

    def return_output(self):
        return self.output_fasta, self.output_csv


class GTF_entry:
    def __init__(self, string):
        self.string = string
        self.values = self.string.split("\t")
        self.id = self.values[0]
        self.start = int(self.values[3])
        self.end = int(self.values[4])
        self.score = float(0.5)  # self.values[5]
        self.sequence = "no sequence set"
        self.length = self.end - self.start

    def __repr__(self):
        return self.sequence[:10] + "..." + f" len={self.length} score={self.score}"

    def set_sequence(self, full_sequence):
        self.sequence = full_sequence[self.start : self.end]


if __name__ == "__main__":
    import argparse

    pass
