def exon_concatenation(
	post_bedtools_fasta: str
) -> list:
	"""Concatenate all sequences starting with identical transcripit ID and outputs it as a list with sequence header (Transcript ID) and concatenated sequences as tuples.

	Args:
		post_bedtools_fasta: The name of the fasta file obtained after bedtools has been run

	Returns:
		A list containing transcript ID and concatenated exons in tuples.
	"""
    with open(post_bedtools_fasta,'r') as fa:
        annotation = []
        fasta_format_list = []
        for line1,line2 in zip(fa,fa):
            if len(annotation) == 0:
                annotation.append(line1[0:16])
                read = line2[:-1]
            else:
                if annotation[-1] == line1[0:16]:
                    read += line2[:-1]
                elif annotation[-1] != line1[0:16]:
                    fasta_format_list.append((annotation[-1],read))
                    annotation.append(line1[0:16])
                    read = line2[:-1]
        fasta_format_list.append((annotation[-1],read))
    return fasta_format_list
