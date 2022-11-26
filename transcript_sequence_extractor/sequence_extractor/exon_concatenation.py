def exon_concatenation_old(
	filename: str
) -> list:
	"""Concatenates all sequences in fasta file with the same transcript ID header and then outputs a list containing sequence headers (Transcript ID) and sequences that have been concatenated.

	Args:
		filename: The name of the fasta file having multiple entries for the same transcript ID.

	Returns:
		A list with headers in the even indices and their corresponding sequences in the odd indices
	"""
	fa = open(filename,'r')
	lines = fa.readlines()
	to_write_to_file = []
	for x in range(int(len(lines)/2)):
		if x == 0:
			annotation = lines[0][0:16]
			read = lines[1][:-1]
		if x >= 1:
			if lines[2*x][1:16] == lines[2*(x-1)][1:16]:
				read+= lines[(2*x)+1][:-1]
			else:
				to_write_to_file.append(annotation)
				to_write_to_file.append(read)
				annotation = lines[2*x][0:16]
				read = lines[(2*x)+1][:-1]
	to_write_to_file.append(annotation)
	to_write_to_file.append(read)
	return to_write_to_file

def exon_concatenation(
	post_bedtools_fasta: str
) -> list:
	"""Concatenate all sequences starting with identical transcripit ID and outputs it as a list with sequence header (Transcript ID) and concatenated sequences.

	Args:
		post_bedtools_fasta: The name of the fasta file obtained after bedtools has been run

	Returns:
		A list with transcript ID in even indices and corresponding concatenated exons in odd indices.
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
                    fasta_format_list.append(annotation[-1])
                    fasta_format_list.append(read)
                    annotation.append(line1[0:16])
                    read = line2[:-1]
        fasta_format_list.append(annotation[-1])
        fasta_format_list.append(read)
    return fasta_format_list
