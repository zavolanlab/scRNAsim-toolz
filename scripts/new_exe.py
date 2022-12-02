import argparse
import transcript_extractor as te
import exon_length_filter as elf
import representative as rtcl
import poisson_sampling as ps
import writegtf as gt
import match_reprtranscript_expressionlevel as ma
def exe(input_file, csv, gtf, input_csv, transcript_nr,input_free = True):
    file_name,source_pathway_name_2,deposit_pathway_name_2 = te.extract_transcript(input_file, Input_free = input_free)
    inter_mediate_file_directory = file_name +"_intermediate_file.txt"
    print("Transcripts are filtered based on transcript score. Please wait...")
    pre_filter_representative_transcripts_dict = rtcl.find_repr_by_SupportLevel(inter_mediate_file_directory)
    print("Transcripts filtered\n")
    dictionary1 = elf.exon_length_filter(file_name,gen_dict= pre_filter_representative_transcripts_dict, Input_free = input_free)
    print(dictionary1)

    tsv_input = ma.match_reprTranscript_expressionLevel(input_csv, dictionary1, inter_mediate_file_directory)
    print("Poisson sampling of transcripts")
    ps.transcript_sampling(transcript_nr, tsv_input, csv)
    print("output csv file ready")
    print("writing output gtf file")
    gt.gtf_file_writer(input_file, csv, gtf)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="transcript sampler",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--annotation", required=True, help="gtf file with genome annotation")
    parser.add_argument("--output_csv", required=True, help="output csv file")
    parser.add_argument("--input_csv", required=True, help="output csv file")
    parser.add_argument("--output_gtf", required=True, help="output gtf file")
    parser.add_argument("--transcript_number", required=True, help="total number of transcripts to sample")
    args = parser.parse_args()
    exe(args.annotation, args.output_csv, args.output_gtf, args.input_csv, args.transcript_number)
