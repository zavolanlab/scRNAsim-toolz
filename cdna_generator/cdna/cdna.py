import sys
import pandas as pd

from gtfparse import read_gtf


# returns GTF with essential columns such as "feature", "seqname", "start", "end"
# alongside the names of any optional keys which appeared in the attribute column
df_input_GTF = read_gtf("Example_GTF_Input.gtf")
df_input_CSV = pd.read_csv("copy_number_input.csv")

df_input_CSV = df_input_CSV.reset_index()  # make sure indexes pair with number of rows

df_input_GTF['Binding_Probability'] = pd.to_numeric(df_input_GTF['Binding_Probability']) # convert to numeric
df_normalization_bind_probablility = df_input_GTF.groupby('seqname')['Binding_Probability'].sum() # extract binding probablility

# Add New columns to the existing DataFrame
df_input_GTF["Normalized_Binding_Probability"] = ''
df_input_GTF["Transcript_Copy_Number"] = ''


# Adds Normalized_Binding_Probability and Transcript_Copy_Number to each transcript in the dataframe
for index, row in df_input_GTF.iterrows():
    # GTF transcript ID 
    id_GTF = str(row['seqname'])                
    # CVS transcript ID 
    id_CSV = str(row['seqname']).split('_')[1]  
    # Calculate Normalized_Binding_Probability and add to GTF dataframe
    df_input_GTF.loc[index, 'Normalized_Binding_Probability'] = row['Binding_Probability'] / df_normalization_bind_probablility[id_GTF]
    # Calculate Normalized_Binding_Probability and add to GTF dataframe
    csv_transcript_copy_number = df_input_CSV.loc[df_input_CSV['ID of transcript'] == int(id_CSV), 'Transcript copy number'].iloc[0]
    df_input_GTF.loc[index,'Transcript_Copy_Number'] = round(csv_transcript_copy_number * df_input_GTF.loc[index,'Normalized_Binding_Probability'])
    

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



if __name__ == "__main__":
    import argparse

    pass
