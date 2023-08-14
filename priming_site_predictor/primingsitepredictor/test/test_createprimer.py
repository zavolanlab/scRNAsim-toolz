"""
Created on Tue Dec 20 14:06:20 2022

@author: RobinC
"""

import os, sys
import pytest
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import createprimer as cp
import unittest


class TestCreatePrimer(unittest.TestCase):

    def test_init(self):
        """Test function for the  __init__() method"""

        # Test default values
        primer = cp.CreatePrimer()
        self.assertEqual(primer.name, 'primer1')
        self.assertEqual(primer.primer_length, 15)
        self.assertEqual(primer.primer_sequence, 'T'*15)
        self.assertEqual(primer.lines, ['<primer1', 'TTTTTTTTTTTTTTT'])

        # Test custom values
        primer = cp.CreatePrimer('my_primer', 20)
        self.assertEqual(primer.name, 'my_primer')
        self.assertEqual(primer.primer_length, 20)
        self.assertEqual(primer.primer_sequence, 'T'*20)
        self.assertEqual(primer.lines, ['<my_primer', 'TTTTTTTTTTTTTTTTTTTT'])


    def test_create_fasta(self):
        """Test function for the  test_create_fasta() method"""

        # Test default values
        primer = cp.CreatePrimer()
        primer.create_fasta()
        with open('primer1.fasta', 'r') as f:
            lines = f.readlines()
        self.assertEqual(lines, ['<primer1\n', 'TTTTTTTTTTTTTTT\n'])

        # Test custom values
        primer = cp.CreatePrimer('my_primer', 20)
        primer.create_fasta()
        with open('my_primer.fasta', 'r') as f:
            lines = f.readlines()
        self.assertEqual(lines, ['<my_primer\n', 'TTTTTTTTTTTTTTTTTTTT\n'])




















