# -*- coding: utf-8 -*-
"""
Changed on Fr Dec 23 14:34:20 2022

@author: RobinC
"""

class CreatePrimer:
    """This Class creates an instance of a primer of a desired length and primer name
    which can be saved as a fasta file. By default the length is 15 and name is primer1"""
    def __init__(self, name='primer1', primerlength=15):
        self.name = name
        self.primer_length = primerlength
        self.primer_sequence = 'T'*self.primer_length
        self.lines = [f'<{self.name}', self.primer_sequence]

    
    def create_fasta(self):
        with open(f'{self.name}.fasta', 'w') as f:
            for line in self.lines:
                f.write(line)
                f.write('\n')




