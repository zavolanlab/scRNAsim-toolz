import random


dna_seq = {
    "ATAACATGTGGATGGCCAGTGGTCGGTTGTTACACGCCTACCGCGATGCTGAATGACCCGGACTAGAGTGGCGAAATTTATGGCGTGTGACCCGTTATGC": 100,
    "TCCATTTCGGTCAGTGGGTCATTGCTAGTAGTCGATTGCATTGCCATTCTCCGAGTGATTTAGCGTGACAGCCGCAGGGAACCCATAAAATGCAATCGTA": 100
}

mean_length = 12
std = 1

term_frags = []
for seq, counts in dna_seq.items():
    for _ in range(counts):
        n_cuts = int(len(seq)/mean_length)
        cuts = random.sample(range(1,len(seq)-1), n_cuts)
        cuts.sort()
        cuts.insert(0,0)
        term_frag = ""
        for i, val in enumerate(cuts):
            if i == len(cuts)-1:
                fragment = seq[val:cuts[-1]]
            else:
                fragment = seq[val:cuts[i+1]]
            if mean_length-std <= len(fragment) <= mean_length+std:
                term_frag = fragment
        if term_frag == "":
            continue
        else:
            term_frags.append(term_frag)
    
with open('terminal_frags.txt', 'w') as f:
    for line in term_frags:
        f.write(line)
        f.write('\n')

