def exon_concatenation:
	fa = open("fasta.fa",'r')
	lines = fa.readlines()
	for x in range(int(len(lines)/2)):
    		if x == 0:
        		annotation = lines[0]
        		read = lines[1]
    		if x >= 1:
        		if lines[2*x] == lines[2*(x-1)]:
            			read+= lines[(2*x)+1]
        		else:
				return annotation
				return read
            			annotation = lines[2*x]
            			read = lines[(2*x)+1]
