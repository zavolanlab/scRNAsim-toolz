def exon_concatenation(filename):
	fa = open(filename,'r')
	lines = fa.readlines()
	to_write_to_file = []
	for x in range(int(len(lines)/2)):
		if x == 0:
			annotation = lines[0]
			read = lines[1]
		if x >= 1:
			if lines[2*x] == lines[2*(x-1)]:
				read+= lines[(2*x)+1]
			else:
				to_write_to_file.append(annotation)
				to_write_to_file.append(read)
				annotation = lines[2*x]
				read = lines[(2*x)+1]
	to_write_to_file.append(annotation)
	to_write_to_file.append(read)
	return to_write_to_file