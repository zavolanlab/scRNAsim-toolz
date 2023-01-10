
def gtf_file_writer (original_file, output_file): 
    output = []
    rep_transcript_dict = get_rep_trans(original_file)

    with open(original_file, 'r') as f:
            for entry in f: 
                if entry[0] != '#':
                    attributes = attributs_converter(entry)
                    type_ = attributes[2]
                    if type_ == 'gene':
                        gene_id = find_in_attributs(attributes, 'gene_id')
                        output.append(entry)
                    if type_ != 'gene':
                        transcript_id = find_in_attributs(attributes, 'transcript_id')
                        if rep_transcript_dict[gene_id] == transcript_id:
                            output.append(entry)

    with open(output_file, 'w') as last_file:
        last_file.write(output)